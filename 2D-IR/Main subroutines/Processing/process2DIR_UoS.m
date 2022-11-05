function dataStruct = process2DIR_UoS(app,dataStruct,ReProcess,bkgIdx,varargin)
% Description:  This function apodises, zeropads and phases 2D IR data
%               from the setup at the University of Sheffield
%
% Usage: dataStruct = process2DIR_UoS(dataStruct)
%        If the manual mode is desired, set debug to 1 and prepare all inputs accordingly
% Data Inputs:
%     cmprobe           (Double)
%     bins              (Double)
%     t2delays          (Double)
%     Ndelays           (Double)
%     Nspectra          (Double)
%     Ndatastates       (Double)
%     Nbins             (Double)
%     Nslowmod          (Double)
%     count             (Cell array)
%     probe             (Cell array)
%     reference         (Cell array)
%     interferogram     (Cell array)
%     signal			(Cell array)
% Additional inputs from the dataStruct structure (GUI options, etc.):
%   ...To be listed
% Outputs:
%   ...To be listed.
% Comments: 
%     It is possible to read both RAW (intensity) and SIGNAL files from the MESS program.
%     The calculation routines were updated accordingly in MESS and now the chopper signal and all
%     datastates should be calculated properly.
%
% Ricardo Fernandez-Teran / 20.04.2022 / v3.0a

%% Choose between debug/manual mode and normal mode
debug=0;
autoplot=0*debug;
show_int=1*debug;

dummy = 'nothing';
% dummy = 'diff'; % Can be the dummy number of 'diff' (will plot 2 - 1) - more options coming soon

%% READ from dataStruct
if isempty(varargin)
    ShowWaitBar = true;
else
    ShowWaitBar = false;
end

%%% Read data
cmprobe         = dataStruct.cmprobe;
est_probe       = dataStruct.est_probe; % Estimated probe from spectrograph info
t1delays        = dataStruct.t1delays;
datatype        = dataStruct.datatype;
Ndelays         = dataStruct.Ndelays;

Ndatastates     = 1;

Nspectra        = dataStruct.Nspectra;
Ndummies        = dataStruct.Ndummies;
Nslowmod        = dataStruct.Nslowmod;
Nbins           = dataStruct.Nbins;
interferogram	= dataStruct.interferogram;
signal		    = dataStruct.signal;
rootdir         = dataStruct.rootdir;
datasubdir      = dataStruct.datafilename;

w0              = dataStruct.w0;
dt1             = dataStruct.dt1;

datadir     = [rootdir filesep datasubdir];

%%% Read GUI options
bkg_sub         = app.I2D_SubtractScatteringCheckBox.Value;
zeropad_enable  = app.I2D_DoZeropaddingCheckBox.Value;
if zeropad_enable == 1
    zeropad_factor  = app.I2D_ZeropadFactor.Value;
else
    zeropad_factor  = 0;
end
zeropad_next2k  = app.I2D_ZeropadToNext2K.Value;
apodise_method  = app.I2D_ApodisationFunction.Value;
probe_calib     = app.I2D_AutocalibrateprobeaxisCheckBox.Value;
pumpcorrection  = app.I2D_PumpcorrectionCheckBox.Value;

%%% Determine whether the data is transient 2D or not, then read the mode
if dataStruct.Transient2D == 1
    Transient2D         = 1;
    Transient2D_mode    = dataStruct.Transient2D_mode;
else
    Transient2D         = 0;
    Transient2D_mode    = 'None';
end

% Hardcoded settings
    DetSz           = dataStruct.DetSz;
    filter_type     = 'Mean';
    filter_points   = 10;
    zeropad_npoints = 1;
    apodize_gaussian= 0;
    apodize_Gcoeff  = 2/3; % Percentage of the Nbins for the decay of the Gaussian
    probe_fitorder  = 2;
% Physical constants
    c_0             = 2.99792458e-5;    % Speed of light in cm/fs

    
% If there is a calibrated WL file in the current EXPDIR, use it
wavenumberfile = 'CalibratedProbe.csv';

if exist([datadir filesep wavenumberfile],'file') == 2
    probe_calib = 2;
    tmp_probe = readmatrix([datadir filesep wavenumberfile]);
elseif exist([rootdir filesep wavenumberfile],'file') == 2
    probe_calib = 2;
    tmp_probe = readmatrix([rootdir filesep wavenumberfile]);
else
    tmp_probe = [];
end

for k=1:length(DetSz)
    if ~isempty(tmp_probe)
        idx = (1:DetSz(k)) + sum(DetSz(1:k-1));
        saved_probe{k}  = tmp_probe(idx);
        cmprobe{k}      = 1:DetSz(k);
    else
        saved_probe{k}  = 1:DetSz(k);
        cmprobe{k}      = 1:DetSz(k);
    end
end 

%% Preprocess the data
% Initialise variables
    FFTinterferogram={}; %#ok<*NASGU> 
    bin=[]; difference={}; coeff={}; scanrange={};
    progress = 0;
    
    
    % Create wait bar if it doesn't exist (when reprocessing data / applying changes)
    if ReProcess == 1 && ShowWaitBar
%         dataStruct.WaitBar = waitbar(0,'Loading data...');
        dataStruct.WBfigure            = uifigure;
        dataStruct.WBfigure.Position(3:4) = [405 175];
        dataStruct.WaitBar             = uiprogressdlg(dataStruct.WBfigure,'Title','2D-IR data processing...','Message','Applying new processing settings...','Icon','info','ShowPercentage','on','Cancelable','on');
    end
    
%% Process the data
for k=1:Nspectra*Ndummies*Nslowmod
for m=1:Ndelays
    if ShowWaitBar
        % Update the Wait Bar
        progress = progress + 1;
    %     waitbar((progress/(Ndelays*Nspectra*Ndummies*2)+0.5),dataStruct.WaitBar,['Processing data... (' num2str(progress) ' of ' num2str(Ndelays*Nspectra*Ndummies) ')']);
        dataStruct.WaitBar.Value    = (progress/(Ndelays*Nspectra*Ndummies*Nslowmod*2)+0.5);
        dataStruct.WaitBar.Message  = ['Processing data... (' num2str(progress) ' of ' num2str(Ndelays*Nspectra*Ndummies*Nslowmod) ')'];
        if dataStruct.WaitBar.CancelRequested
            delete(dataStruct.WBfigure);
            error('User aborted loading the data!');
        end
        drawnow;
    end
    % Remove background (and slow fluctuations in the data, if selected)
    if ReProcess == 0
        switch filter_type
            case 'Median'
                interferogram{m,k}  = -(interferogram{m,k}-mean(interferogram{m,k}));
                interferogram{m,k}  = interferogram{m,k}-medfilt1(interferogram{m,k},filter_points);
                signal{m,k}         = signal{m,k}-mean(signal{m,k},1);
                signal{m,k}         = signal{m,k}-medfilt1(signal{m,k},filter_points);
            case 'Mean'
                signal{m,k}         = signal{m,k}-mean(signal{m,k},1);
        end
    end

binzero{m,k}    = 1; % Because we are using the shaper
binspecmax(m,k) = 1;    
%% Apodise and phase the data
% Apodise the data by the selected method
    % Cosine
    switch apodise_method
        case '0'
            cos_exp=0;
        case '1'
            cos_exp=1;
        case '2'
            cos_exp=2;
        case '3'
            cos_exp=3;
    end
    q=1:Nbins;
    cosine{m,k}             = ((cos(pi*(q-binzero{m,k})/(2*(Nbins-binzero{m,k})))).^cos_exp)';
    box{m,k}                = (heavisideRF(q-binzero{m,k}))';

    % Gaussian
    if apodize_gaussian==1
        gaussian{m,k}           = (exp(-((q-binzero{m,k})./(apodize_Gcoeff.*Nbins)).^2))';
        apodize_function{m,k}   = cosine{m,k}.*gaussian{m,k};
    else
        apodize_function{m,k}   = cosine{m,k};
    end

% Multiply the interferogram by the apodization function (No box! - Use both sides)
    apo_interferogram{m,k}      = interferogram{m,k};
    
% Multiply the signal by the apodization function (Use box! - discard data before BinZero + 1st point by 1/2)
    apo_signal{m,k}             = signal{m,k}.*apodize_function{m,k}.*box{m,k};
    apodize_function{m,k}       = apodize_function{m,k}.*box{m,k};
    
% Decide the FT size
    N_FTpoints{m,k}             = zeropad_factor*Nbins;
    if zeropad_next2k==1
        N_FTpoints{m,k}         = 2^nextpow2(N_FTpoints{m,k});
    end
    
% Calculate the pump frequency axis from interferogram
    Resolution(m,k)             = 1/(N_FTpoints{m,k}*dt1*c_0);
    PumpAxis{m,k}               = ((1:1:N_FTpoints{m,k})-1)'.*Resolution(m,k)+w0;
    
%% Zeropad the data
% Get the points for zeropadding
if apodise_method=="Box"
    mean_interf{m,k}= mean(interferogram{m,k}(end-zeropad_npoints:end),1);
    mean_sig{m,k}   = mean(signal{m,k}(end-zeropad_npoints:end,:),1);
else
    mean_interf{m,k}= 0;
    mean_sig{m,k}   = zeros(1,length(cmprobe{k}));
end

% Build the zeropad matrices
    N_newpoints{m,k}            = N_FTpoints{m,k}-Nbins;
    pad_interf{m,k}             = ones(N_newpoints{m,k},1)*mean_interf{m,k};
    pad_signal{m,k}             = ones(N_newpoints{m,k},length(cmprobe{k})).*mean_sig{m,k};

% Do the zeropadding
    zeropad_interf{m,k}         = [apo_interferogram{m,k};pad_interf{m,k}];
    zeropad_signal{m,k}         = [apo_signal{m,k};pad_signal{m,k}];
        
%% Calculate and fit the phase        
% Shift the interferogram to BinZero and get the phase
%     phased_ZPint{m,k}           = circshift(zeropad_interf{m,k},-(binzero{m,k}));
%     phased_ZPsig{m,k}           = circshift(zeropad_signal{m,k},-(binzero{m,k}),1);

    phased_ZPint{m,k}           = zeropad_interf{m,k};
    phased_ZPsig{m,k}           = zeropad_signal{m,k};
    
    absFFT_ZPint{m,k}           = phased_ZPint{m,k};
    
% Calculate the FFT and phase
    FFT_ZPint{m,k}              = phased_ZPint{m,k};
    FFT_ZPsig{m,k}              = fft(phased_ZPsig{m,k});
    
    ZP_phase{m,k}               = zeros(size(phased_ZPint{m,k}));
    
% The phase is a flat zero because of th shaper
    fittedPhase{m,k}    = ZP_phase{m,k};
    phase_coeff{m,k}    = NaN;

% Phase the data (MCT data + interferogram)
    phasingterm{m,k}            = exp(-1i.*fittedPhase{m,k});
    phased_FFTZPint{m,k}        = FFT_ZPint{m,k}.*phasingterm{m,k};

% Do pump correction
% if pumpcorrection == 1
%     phased_FFTZPsig{m,k}        = real(FFT_ZPsig{m,k}.*phasingterm{m,k})./(abs(phased_FFTZPint{m,k}*ones(1,length(cmprobe{k}))));
%     magnitude_FFTsig{m,k}       = abs(FFT_ZPsig{m,k}-mean(FFT_ZPsig{m,k}))./(real(phased_FFTZPint{m,k}*ones(1,length(cmprobe{k}))));
% else
    phased_FFTZPsig{m,k}        = real(FFT_ZPsig{m,k}.*phasingterm{m,k});
    magnitude_FFTsig{m,k}       = abs(FFT_ZPsig{m,k});
% end
    
%% PROBE AXIS CALIBRATION USING SCATTERING - LINEAR FIT OF THE MAXIMA ALONG THE DIAGONAL
if m==bkgIdx
      % Get a list of the maxima of the interferogram (N_FTpoints/50 around the spectral max) 
      %%%! TO DO: needs to be redefined in terms of cm-1
      switch probe_calib 
            % 0=No calibration, use generated probe axis
            % 1=Doing calibration from scattering, use the results (TBD)
            % 2=Using saved probe (in file)
            case {0,1}
                if unique(diff(est_probe{k})) == 1
                    probe    = linspace(min(PumpAxis{m,k}),max(PumpAxis{m,k}),DetSz(k));
                else
                    probe    = est_probe{k};
                end
            case 2
                probe        = saved_probe{k};
      end

        minProbe_est_idx        = findClosestId2Val(PumpAxis{m,k},min(probe));  
        maxProbe_est_idx        = findClosestId2Val(PumpAxis{m,k},max(probe));          
        fitrange                = minProbe_est_idx:maxProbe_est_idx;
        Npixels                 = size(magnitude_FFTsig{m,k},2);
        pixels                  = transpose(1:1:Npixels);
        % Get the maxima for probe calibration    
        magFFT                  = abs(FFT_ZPsig{m,k});
        [~,maxindex]            = max(magFFT(fitrange,:));
        scattering_maxima_d     = PumpAxis{m,k}(fitrange(maxindex));
%        % Do the fit [OLD WAY]
%         freq_coeff              = polyfit(pixels,scattering_maxima,probe_fitorder);
%         freq_fit                = transpose(polyval(freq_coeff,pixels));
        
        % Do the fit [NEW WAY - ROBUST]
        switch probe_fitorder
            case 1
                model       = 'poly1';
            case 2
                model       = 'poly2';
        end
        warning('off','curvefit:fit:iterationLimitReached');
        mdl                     = fit(pixels,scattering_maxima_d,model,'Robust','Bisquare');
        freq_fit_tmp            = mdl(pixels);
%       % Decide which kind of probe axis to use. It will anyway store all of them
%         switch probe_calib
%             case 1
%             % Check if the Probe Axis makes sense, otherwise use either the stored one
%                 if issorted(freq_fit)
%                     ProbeAxis{k}  = freq_fit;
%                 else
%                     ProbeAxis{k}  = cmprobe{k};
%                 end
%             case 0
%                 ProbeAxis{k}      = cmprobe{k};
%             case 2
%                 ProbeAxis{k}      = saved_probe{k};
%         end

freq_fit{k}         = freq_fit_tmp;   %#ok<*AGROW> 
scattering_maxima{k}= scattering_maxima_d;

switch probe_calib 
    % 0=No calibration, use generated probe axis
    % 1=Doing calibration from scattering, use the results (TBD)
    % 2=Using saved probe (in file)
    case 1
    % Check if the Probe Axis makes sense, otherwise use either the stored one
        if issorted(freq_fit_tmp)
            ProbeAxis{k}    = freq_fit_tmp;
        else
            ProbeAxis{k}    = est_probe{k};
        end
    case 0
        ProbeAxis{k}        = est_probe{k};
    case 2
        ProbeAxis{k}        = saved_probe{k};
end

end
%% Subtract the scattering background (if enabled) and correct the sign of the signals
% If DEBUG is OFF, don't consider the sign of the pump  
if debug==0
    SignPump(m,k)=1;
elseif debug==1 % If debug is ON, consider the sign of the pump
    SignPump(m,k)=sign(real(phased_FFTZPint{m,k}(binspecmax(m,k))));
end

end % Ndelays
end % Ndatastates

for k=1:Nspectra*Ndummies*Nslowmod
for m=1:Ndelays
    % Background subtraction
    if bkg_sub==1 && m~=bkgIdx
        PROC_2D_DATA{m,k}       = (phased_FFTZPsig{m,k}-phased_FFTZPsig{bkgIdx,k})*SignPump(m,k);
    else
        PROC_2D_DATA{m,k}       = phased_FFTZPsig{m,k}*SignPump(m,k);
    end
end
end
% %% SPECTROMETER CALIBRATION
% 
% % Get the maxima for probe calibration    
% magFFT                  = abs(FFT_ZPsig{m,k});
% [~,maxindex]            = max(magFFT(fitrange,:));
% scattering_maxima       = PumpAxis{m,k}(fitrange(maxindex));
% 
% % Do the fit
% switch probe_fitorder
%     case 1
%         model       = 'poly1';
%     case 2
%         model       = 'poly2';
% end
% warning('off','curvefit:fit:iterationLimitReached');
% mdl                     = fit(pixels,scattering_maxima,model,'Robust','Bisquare');
% freq_fit                = mdl(pixels);
% % Decide which kind of probe axis to use. It will anyway store all of them
% switch probe_calib
%     case 1
%     % Check if the Probe Axis makes sense, otherwise use either the stored one
%         if issorted(freq_fit)
%             ProbeAxis       = freq_fit;
%         else
%             ProbeAxis       = cmprobe;
%         end
%     case 0
%         ProbeAxis       = cmprobe;
%     case 2
%         ProbeAxis       = saved_probe;
% end
%     
% %% Transient 2D IR processing
% % % Determine whether it is a transient 2D dataset or not, then do the stuff
% if Transient2D
%     switch Transient2D_mode
%         case 'Spectrum 1 (UV on)' % Show spectrum 1
%             PROCDATA = PROC_2D_DATA(:,1);
%         case 'Spectrum 2 (UV off)' % Show spectrum 2
%             PROCDATA = PROC_2D_DATA(:,2);
%         case 'Difference (1-2)' % Show difference 1-2
%             for m=1:Ndelays
%                 PROCDATA{m,1} = PROC_2D_DATA{m,1} - PROC_2D_DATA{m,2};
%             end
%         case 'Difference (2-1)'
%             for m=1:Ndelays
%                 PROCDATA{m,1} = PROC_2D_DATA{m,2} - PROC_2D_DATA{m,1};
%             end
%         case 'Sum' % Show sum
%             for m=1:Ndelays
%                 PROCDATA{m,1} = PROC_2D_DATA{m,1} + PROC_2D_DATA{m,2};
%             end
%         otherwise
%             return
%     end
%     PROC_2D_DATA    = PROCDATA;
% end

%% WRITE to dataStruct
    dataStruct.ProbeAxis           = ProbeAxis;
    dataStruct.freq_fit            = [];
    dataStruct.PumpAxis            = PumpAxis;
    
    dataStruct.scattering_maxima   = scattering_maxima;
    dataStruct.freq_fit            = freq_fit;
    
    dataStruct.t1delays            = t1delays;
    dataStruct.interferogram       = interferogram;
    dataStruct.binzero             = binzero;
    dataStruct.binspecmax          = binspecmax;
    dataStruct.apodize_function    = apodize_function;
    dataStruct.FFT_ZPsig           = FFT_ZPsig;
    dataStruct.phased_FFTZPsig     = phased_FFTZPsig;
	dataStruct.phased_FFTZPint     = phased_FFTZPint;
    dataStruct.fittedPhase         = fittedPhase;
    dataStruct.phasepoints         = [];
    dataStruct.ZP_phase            = ZP_phase;
    dataStruct.phase_coeff         = phase_coeff;
    dataStruct.apo_interferogram   = apo_interferogram;
    dataStruct.apo_signal          = apo_signal;
    dataStruct.signal              = signal;
    dataStruct.PROC_2D_DATA        = PROC_2D_DATA;
    dataStruct.SpecDiff            = 0;

%%%%%%%% EOF