function dataStruct = process2DIR(app,dataStruct,ReProcess)
% Description: This function apodizes, zeropads and phases 2D IR data
%
% Usage: dataStruct = process2DIR(dataStruct)
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
% Ricardo Fernández-Terán / 22.08.2018 / v3.5a

%% Choose between debug/manual mode and normal mode
debug=0;
autoplot=0*debug; show_int=1*debug;

dummy = 1; % Can be the dummy number of 'diff' (will plot 2 - 1) - more options coming soon

%% READ from dataStruct

% Read data
cmprobe         = dataStruct.cmprobe;
t1delays        = dataStruct.t1delays;
datatype        = dataStruct.datatype;
Ndelays         = dataStruct.Ndelays;
if strcmp(datatype,'Raw')
    Ndatastates = dataStruct.Ndatastates;
else
    Ndatastates = 1;
end
Nspectra        = dataStruct.Nspectra;
Ndummies        = dataStruct.Ndummies;
Nbins           = dataStruct.Nbins;
interferogram	= dataStruct.interferogram;
signal		    = dataStruct.signal;
rootdir         = dataStruct.rootdir;

% Read GUI options
bkg_sub         = app.I2D_SubtractScatteringCheckBox.Value;
zeropad_enable  = app.I2D_DoZeropaddingCheckBox.Value;
if zeropad_enable == 1
    zeropad_factor  = app.I2D_ZeropadFactor.Value;
else
    zeropad_factor  = 0;
end
% Determine whether the data is transient 2D or not, then read the mode
if dataStruct.Transient2D == 1
    Transient2D         = 1;
    Transient2D_mode    = dataStruct.Transient2D_mode;
else
    Transient2D         = 0;
    Transient2D_mode    = 0;
end

zeropad_next2k  = app.I2D_ZeropadToNext2K.Value;
apodise_method  = app.I2D_ApodisationFunction.Value;
phase_method    = app.I2D_PhaseFitFunction.Value;
phase_points    = app.I2D_PhaseFitRangeEdit.Value;
probe_calib     = app.I2D_AutocalibrateprobeaxisCheckBox.Value;
binzero         = num2cell(-1*ones(Ndelays,Ndatastates*Nspectra*Ndummies));
pumpcorrection  = app.I2D_PumpcorrectionCheckBox.Value;
% If there is a calibrated WL file in the current ROOTDIR, use it
if probe_calib == 0 && exist([rootdir filesep 'CalibratedProbe.csv'],'file') == 2
    probe_calib = 2;
    wavenumberfile='CalibratedProbe.csv';
    saved_probe=csvread([rootdir filesep wavenumberfile]);
end    

% Hardcoded settings
    filter_type     = 'Median';
    filter_points   = 10;
    zeropad_npoints = 10;
    apodize_gaussian= 0;
    apodize_Gcoeff  = 2/3; % Percentage of the Nbins for the decay of the Gaussian
    probe_fitorder  = 2;
    phase_fitmethod = 'Shift bins'; % 'Shift wavenumbers' or 'Shift bins'
    phase_wavenum   = 100;
    phase_shiftbins = 10;
% Physical constants
    HeNe            = 2.11079;          % HeNe period (fs)
    c_0             = 2.99792458e-5;    % Speed of light in cm/fs
    
%% Preprocess the data
% Initialise variables
    FFTinterferogram={}; binspecmax=[]; bininterfmax=[];
    bin=[]; difference={}; coeff={}; scanrange={};
    progress = 0;
    
    % Create wait bar if it doesn't exist (when reprocessing data / applying changes)
    if ReProcess == 1
%         dataStruct.WaitBar = waitbar(0,'Loading data...');
        dataStruct.WBfigure            = uifigure;
        dataStruct.WBfigure.Position(3:4) = [405 175];
        dataStruct.WaitBar             = uiprogressdlg(dataStruct.WBfigure,'Title','2D-IR data processing...','Message','Loading data...','Icon','info','ShowPercentage','on','Cancelable','on');
    end
    
%% Process the data
for k=1:Nspectra*Ndummies
for m=1:Ndelays
    % Update the Wait Bar
    progress = progress + 1;
%     waitbar((progress/(Ndelays*Nspectra*Ndummies*2)+0.5),dataStruct.WaitBar,['Processing data... (' num2str(progress) ' of ' num2str(Ndelays*Nspectra*Ndummies) ')']);
    dataStruct.WaitBar.Value    = (progress/(Ndelays*Nspectra*Ndummies*2)+0.5);
    dataStruct.WaitBar.Message  = ['Processing data... (' num2str(progress) ' of ' num2str(Ndelays*Nspectra*Ndummies) ')'];
    if dataStruct.WaitBar.CancelRequested
        delete(dataStruct.WBfigure);
        error('User aborted loading the data!');
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
                interferogram{m,k}  = -(interferogram{m,k}-mean(interferogram{m,k}));
                signal{m,k}         = signal{m,k}-mean(signal{m,k},1);
        end
    end
    FFTinterferogram{m,k}       = fft(interferogram{m,k});
    Nbins                       = length(interferogram{m,k});

%% Find BinZero before apodisation
% Find BinSpecMax and BinInterfMax
    [~,bininterfmax(m,k)]       = max(interferogram{m,k});
    [~,binspecmax(m,k)]         = max(abs(FFTinterferogram{m,k}(20:floor(length(FFTinterferogram{m,k})/2))));
    binspecmax(m,k)             = binspecmax(m,k) + 20;
    Resolution(m,k)             = 1/(Nbins*c_0*HeNe);
    switch phase_fitmethod
        case 'Shift wavenumbers'
            shift_findphase(m,k)        = round(phase_wavenum/Resolution(m,k));
        case 'Shift bins'
            shift_findphase(m,k)        = phase_shiftbins;
        otherwise
            shift_findphase(m,k)        = 10;
    end
% If binzero is not known (= -1) then try to guess it
    if binzero{m,k}==-1
        % The method will try shifting the interferogram until the point where the phase is flat
        for p=1:200
            bin{m,k}(p)         = p+bininterfmax(m,k)-100;
            % Unwrap the phases to remove discontinuities greater than 2*pi
            tempPhase{m,k}      = unwrap(angle(fft(circshift(interferogram{m,k},-bin{m,k}(p)))));
            difference{m,k}(p)  = (tempPhase{m,k}(binspecmax(m,k)+shift_findphase(m,k)))-(tempPhase{m,k}(binspecmax(m,k)-shift_findphase(m,k)));
        end
        % Then, it will fit a straight line that tells how much the slope is changing
        coeff{m,k}              = polyfit(bin{m,k},difference{m,k},1);
        % Afterwards, it will find the zero crossing point from the fitted coefficients
        binzero{m,k}            = round(-(coeff{m,k}(2)/coeff{m,k}(1)));
    end
    
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
    box{m,k}                = (heaviside(q-binzero{m,k}))';

    % Gaussian
    if apodize_gaussian==1
        gaussian{m,k}           = (exp(-((q-binzero{m,k})./(apodize_Gcoeff.*Nbins)).^2))';
        apodize_function{m,k}   = cosine{m,k}.*gaussian{m,k};
    else
        apodize_function{m,k}   = cosine{m,k};
    end

% Multiply the interferogram by the apodization function (No box! - Use both sides)
    apo_interferogram{m,k}      = interferogram{m,k}.*apodize_function{m,k};
    
% Multiply the signal by the apodization function (Use box! - discard data before BinZero + 1st point by 1/2)
    apo_signal{m,k}             = signal{m,k}.*apodize_function{m,k}.*box{m,k};
    apodize_function{m,k}       = apodize_function{m,k}.*box{m,k};
    
% Decide the FT size
    N_FTpoints{m,k}             = zeropad_factor*Nbins;
    if zeropad_next2k==1
        N_FTpoints{m,k}         = 2^nextpow2(N_FTpoints{m,k});
    end
    
% Calculate the pump frequency axis from interferogram by using HeNe counting
    Resolution(m,k)             = 1/(N_FTpoints{m,k}*HeNe*c_0);
    PumpAxis{m,k}               = ((1:1:N_FTpoints{m,k})-1)'.*Resolution(m,k);
    
%% Zeropad the data
% Get the points for zeropadding
if apodise_method=="Box"
    mean_interf{m,k}= mean(interferogram{m,k}(end-zeropad_npoints:end),1);
    mean_sig{m,k}   = mean(signal{m,k}(end-zeropad_npoints:end,:),1);
else
    mean_interf{m,k}= 0;
    mean_sig{m,k}   = zeros(1,length(cmprobe));
end

% Build the zeropad matrices
    N_newpoints{m,k}            = N_FTpoints{m,k}-Nbins;
    pad_interf{m,k}             = ones(N_newpoints{m,k},1)*mean_interf{m,k};
    pad_signal{m,k}             = ones(N_newpoints{m,k},length(cmprobe)).*mean_sig{m,k};

% Do the zeropadding
    zeropad_interf{m,k}         = [apo_interferogram{m,k};pad_interf{m,k}];
    zeropad_signal{m,k}         = [apo_signal{m,k};pad_signal{m,k}];
        
%% Calculate and fit the phase        
% Shift the interferogram to BinZero and get the phase
    phased_ZPint{m,k}           = circshift(zeropad_interf{m,k},-(binzero{m,k}));
    phased_ZPsig{m,k}           = circshift(zeropad_signal{m,k},-(binzero{m,k}),1);
       
    absFFT_ZPint{m,k}           = abs(fft(phased_ZPint{m,k}));
    [~,binspecmax(m,k)]         = max(absFFT_ZPint{m,k}(100:floor(length(absFFT_ZPint{m,k})/2)));
    binspecmax(m,k)             = binspecmax(m,k) + 100;
    
% Calculate the FFT and phase
    FFT_ZPint{m,k}              = fft(phased_ZPint{m,k});
    FFT_ZPsig{m,k}              = fft(phased_ZPsig{m,k});
    
    ZP_phase{m,k}               = unwrap(angle(FFT_ZPint{m,k}));
    
% Fit the phase to the whole interferogram
    shift_fitphase(m,k)         = round(phase_points/Resolution(m,k));
    points{m,k}                 = transpose((binspecmax(m,k)-shift_fitphase(m,k)):1:(binspecmax(m,k)+shift_fitphase(m,k)));
    
    switch phase_method
        case 'Constant'
            fittedPhase{m,k}    = ones(length(ZP_phase{m,k}),1)*mean(ZP_phase{m,k}(points{m,k}));
            phase_coeff{m,k}    = mean(ZP_phase{m,k}(points{m,k}));
        case 'Linear'
            phase_coeff{m,k}    = polyfit(points{m,k},ZP_phase{m,k}(points{m,k}),1);
            fittedPhase{m,k}    = transpose(polyval(phase_coeff{m,k},1:1:N_FTpoints{m,k}));
        case 'Quadratic'
            phase_coeff{m,k}    = polyfit(points{m,k},ZP_phase{m,k}(points{m,k}),2);
            fittedPhase{m,k}    = transpose(polyval(phase_coeff{m,k},1:1:N_FTpoints{m,k}));
        case 'Cubic'
            phase_coeff{m,k}    = polyfit(points{m,k},ZP_phase{m,k}(points{m,k}),3);
            fittedPhase{m,k}    = transpose(polyval(phase_coeff{m,k},1:1:N_FTpoints{m,k}));
        case 'No fit'
            fittedPhase{m,k}    = ZP_phase{m,k};
            phase_coeff{m,k}    = NaN;
    end

% Phase the data (MCT data + interferogram)
    phasingterm{m,k}            = exp(-1i*fittedPhase{m,k});
    phased_FFTZPint{m,k}        = FFT_ZPint{m,k}.*phasingterm{m,k};

% Do pump correction
if pumpcorrection == 1
    phased_FFTZPsig{m,k}        = real(FFT_ZPsig{m,k}.*phasingterm{m,k})./(abs(phased_FFTZPint{m,k}*ones(1,length(cmprobe))));
    magnitude_FFTsig{m,k}       = abs(FFT_ZPsig{m,k}-mean(FFT_ZPsig{m,k}))./(real(phased_FFTZPint{m,k}*ones(1,length(cmprobe))));
else
    phased_FFTZPsig{m,k}        = real(FFT_ZPsig{m,k}.*phasingterm{m,k});
    magnitude_FFTsig{m,k}       = abs(FFT_ZPsig{m,k});
end
    
%% PROBE AXIS CALIBRATION USING SCATTERING - LINEAR FIT OF THE MAXIMA ALONG THE DIAGONAL
if m==1 && k==1
      % Get a list of the maxima of the interferogram (N_FTpoints/50 around the spectral max) 
      %%%! TO DO: needs to be redefined in terms of cm-1
        P                       = floor(N_FTpoints{m,k}/50);
        fitrange                = (binspecmax(m,k)-P):(binspecmax(m,k)+P);
        Npixels                 = size(magnitude_FFTsig{m,k},2);
        pixels                  = transpose(1:1:Npixels);
        % Get the maxima for probe calibration    
        magFFT                  = abs(FFT_ZPsig{m,k});
        [~,maxindex]            = max(magFFT(fitrange,:));
        scattering_maxima       = PumpAxis{m,k}(fitrange(maxindex));
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
        mdl                     = fit(pixels,scattering_maxima,model,'Robust','Bisquare');
        freq_fit                = mdl(pixels);
      % Decide which kind of probe axis to use. It will anyway store all of them
        switch probe_calib
            case 1
            % Check if the Probe Axis makes sense, otherwise use either the stored one
                if issorted(freq_fit)
                    ProbeAxis       = freq_fit;
                else
                    ProbeAxis       = cmprobe;
                end
            case 0
                ProbeAxis       = cmprobe;
            case 2
                ProbeAxis       = saved_probe;
        end
      
end
%% Subtract the scattering background (if enabled) and correct the sign of the signals
% If DEBUG is OFF, don't consider the sign of the pump  
if debug==0
    SignPump(m,k)=1;
elseif debug==1 % If debug is ON, consider the sign of the pump
    SignPump(m,k)=sign(real(phased_FFTZPint{m,k}(binspecmax(m,k))));
end

% Now do the stuff
if bkg_sub==1 && m~=1
    PROC_2D_DATA{m,k}       = (phased_FFTZPsig{m,k}-phased_FFTZPsig{1,1})*SignPump(m,k);
else
    PROC_2D_DATA{m,k}       = phased_FFTZPsig{m,k}*SignPump(m,k);
end

end % Ndelays
end % Ndatastates

%% Transient 2D IR processing
% % Determine whether it is a transient 2D dataset or not, then do the stuff
if Transient2D
    switch Transient2D_mode
        case 1 % Show spectrum 1
            PROCDATA = PROC_2D_DATA(:,1);
        case 2 % Show spectrum 2
            PROCDATA = PROC_2D_DATA(:,2);
        case 3 % Show difference 1-2
            for m=1:Ndelays
                PROCDATA{m,1} = PROC_2D_DATA{m,1} - PROC_2D_DATA{m,2};
            end
        case 4 % Show difference 2-1
            for m=1:Ndelays
                PROCDATA{m,1} = PROC_2D_DATA{m,2} - PROC_2D_DATA{m,1};
            end
        case 5 % Show sum
            for m=1:Ndelays
                PROCDATA{m,1} = PROC_2D_DATA{m,1} + PROC_2D_DATA{m,2};
            end
    end
    PROC_2D_DATA    = PROCDATA;
end

%% Dummies processing
if isnumeric(dummy)
    PROC_2D_DATA        = PROC_2D_DATA(:,dummy);
%     warndlg(['Plotting dummy ' num2str(dummy)]);
elseif strcmp(dummy,'diff')
    for i=1:Ndelays
        PROCDATA{i,1}   = PROC_2D_DATA{i,2} - PROC_2D_DATA{i,1};
    end
    PROC_2D_DATA        = PROCDATA;
%     warndlg('Plotting dummy 2 - dummy 1');
end

% interferogram   = 
%% WRITE to dataStruct
    dataStruct.ProbeAxis           = ProbeAxis;
    dataStruct.freq_fit            = freq_fit;
    dataStruct.scattering_maxima   = scattering_maxima;
    dataStruct.PumpAxis            = PumpAxis;
    
    dataStruct.t1delays            = t1delays;
    dataStruct.interferogram       = interferogram;
    dataStruct.binzero             = binzero;
    dataStruct.binspecmax          = binspecmax;
    dataStruct.apodize_function    = apodize_function;
    dataStruct.FFT_ZPsig           = FFT_ZPsig;
    dataStruct.phased_FFTZPsig     = phased_FFTZPsig;
	dataStruct.phased_FFTZPint		= phased_FFTZPint;
    dataStruct.fittedPhase			= fittedPhase;
    dataStruct.phasepoints         = points;
    dataStruct.ZP_phase            = ZP_phase;
    dataStruct.phase_coeff         = phase_coeff;
    dataStruct.apo_interferogram   = apo_interferogram;
    dataStruct.apo_signal          = apo_signal;
    dataStruct.signal              = signal;
    dataStruct.PROC_2D_DATA        = PROC_2D_DATA;
    dataStruct.SpecDiff            = 0;
%% Stuff to do only if DEBUG MODE is ON
if debug==1
    % Export the data properly
        dataStruct.PROC_2D_DATA{1,1}   = PROC_2D_DATA{1,1};
%         dataStruct.t2_delay_ps         = t2_delay_ps;
        dataStruct.t2_delay_ps         = 0;
    % Show phased interferogram
    if show_int ==1
        figure
        plot(PumpAxis{1,1},real(phased_FFTZPint{1,1}*SignPump(1,1)))
        hl=refline(0,0); hl.Color = [0.5 0.5 0.5];
        xlim([min(ProbeAxis)*0.75,max(ProbeAxis)*1.25]);
    end
    % Make 2D plot
    if autoplot==1
        figure
        plotaxis = gca; plot_2DIR(dataStruct,plotaxis);
    end
end