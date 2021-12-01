function  dataStruct = load2DIRlab3(dataStruct,CurrDelay)
PowerSpectrum=0;
% Description: This function loads all 2DIR data in the file format from the Lab 3 program
% Usage: dataStruct = load2DIRlab3(dataStruct)
% Inputs:
%   dataStruct structure with fields:
%     datafilename
%     rootdir
%     
% Outputs:
%     cmprobe           (Double)
%     bins              (Double)
%     t2delays          (Double)
%     Ndelays           (Double)
%     Nspectra          (Double)
%     Ndatastates       (Double)
%     Nbins             (Double)
%     Nslowmod          (Double)
%     interferogram     (Cell array)
%     signal			(Cell array)
%
% Ricardo Fernández-Terán / 17.12.2018 / v1.0c

% Read options from GUI

% Read from dataStruct
Ndelays             = dataStruct.Ndelays;
rootdir             = char(dataStruct.rootdir);
filelist            = dir(rootdir);
filenames           = {filelist.name}';
filenames           = filenames(~[filelist.isdir]);
if CurrDelay == 0
    filename        = [rootdir filesep filenames{contains(filenames,'Pinhole','IgnoreCase',true)}];
else
    filenames_delays    = filenames(~contains(filenames,'Pinhole','IgnoreCase',true));
    filenames_delays    = filenames_delays(~contains(filenames_delays,'Phases','IgnoreCase',true));
    t2delays            = dataStruct.t2delays;
    [dataStruct.t2delays,t2Index]    = sort(t2delays);
    filenames_delays    = filenames_delays(t2Index);
    filename            = filenames_delays{CurrDelay};
    filename            = [rootdir filesep filename];
end

phasefile       = [rootdir filesep 'Phases.txt'];
zeropad_next2k  = dataStruct.zeropad_next2k.Value;
apodise_method  = char(dataStruct.apodise_method.String{dataStruct.apodise_method.Value});
pumpcorrection  = dataStruct.PumpCorrection_tick.Value;

% Read GUI options
bkg_sub         = dataStruct.BkgSubTick.Value;
zeropad_enable  = dataStruct.zeropad_tick.Value;

% Hardcoded settings
    filter_type     = 'Median';
    filter_points   = 10;
    zeropad_npoints = 10;
    probe_fitorder  = 2;

% Physical constants
    HeNe        = 2.11079;          % HeNe period (fs)
    c_0         = 2.99792458e-5;    % Speed of light in cm/fs

% Open the phasefile (if it exists)
if exist(phasefile,'file') ~= 0
    phases  = dlmread(phasefile,'\t');
else
    phases  = pi*ones(Ndelays);
end

% Open the files for the corresponding T2 delay
ds1         = importdata(filename,',',1);
ds2         = importdata(filename,',',36);
ds3         = importdata(filename,',',71);
ds4         = importdata(filename,',',106);

% Join the four states in a 3D array
data(:,:,1) = ds1.data(1:end-1,:)';
data(:,:,2) = ds2.data(1:end-1,:)';
data(:,:,3) = ds3.data(1:end-1,:)';
data(:,:,4) = ds4.data(1:end-1,:)';

% Initialize some variables
Nbins       = 0.5*(size(data,1)+1);
Npixels     = size(data,2);
t1delays    = linspace(0,Nbins-1,Nbins)*HeNe;

% Decide about zeropadding
if zeropad_enable == 1
    zeropad_factor  = str2double(dataStruct.zeropad_factor.String);
    if zeropad_next2k
        zeropad_factor = 2^nextpow2(Nbins*zeropad_factor)/Nbins;
    end
else
    zeropad_factor  = 1;
end

% Split the data into rephasing and non-rephasing
data_NR     = data(1:Nbins,:,:);
data_RE     = data(Nbins:end,:,:);

% Add the datastates to obtain the correct signal
% CORRECT ADDITION: 1 - 2 + 3 - 4
% signal_xx = Nbins x Npixels
signal_NR   = data_NR(:,:,1) - data_NR(:,:,2) + data_NR(:,:,3) - data_NR(:,:,4);
signal_RE   = data_RE(:,:,1) - data_RE(:,:,2) + data_RE(:,:,3) - data_RE(:,:,4);

% Calculate W1 axis
PumpAxis = ((1:1:Nbins*zeropad_factor)-1)./(Nbins*zeropad_factor*HeNe*c_0);

% Show wait bar
waitbar(CurrDelay/Ndelays,dataStruct.WaitBar,['Processing data... (' num2str(CurrDelay) ' of ' num2str(Ndelays) ')']);
    

if CurrDelay == 0 % This is the pinhole. Use it to calibrate the W3 axis
    Power_spec              = abs(fft(signal_NR,Nbins*zeropad_factor))+abs(fft(signal_RE,Nbins*zeropad_factor));
    [~,binspecmax]          = max(sum(Power_spec,2));
    P                       = floor(Nbins/50);
    fitrange                = (binspecmax-P):(binspecmax+P);
    pixels                  = transpose(1:1:Npixels);
    % Get the maxima for probe calibration    
    [~,maxindex]            = max(Power_spec,[],1);
    scattering_maxima       = PumpAxis(maxindex)';
    switch probe_fitorder
        case 1
            model       = 'poly1';
        case 2
            model       = 'poly2';
    end
    warning('off','curvefit:fit:iterationLimitReached');
    mdl                     = fit(pixels,scattering_maxima,model,'Robust','Bisquare');
    ProbeAxis               = mdl(pixels);
    %% WRITE TO dataStruct
    dataStruct.freq_fit            = ProbeAxis;
    dataStruct.scattering_maxima   = scattering_maxima;
    dataStruct.interferogram       = signal_RE(:,16); %
    dataStruct.binspecmax_pinhole  = binspecmax;%
    dataStruct.ProbeAxis           = ProbeAxis;
    dataStruct.cmprobe             = ProbeAxis;
    dataStruct.powerspec           = Power_spec;
else    % This is a sample t2 delay, process it normally
    % Apodize the data
    switch apodise_method
        case 'Box'
            cos_exp=0;
        case 'Cos'
            cos_exp=1;
        case 'Cos^2'
            cos_exp=2;
        case 'Cos^3'
            cos_exp=3;
    end
    q=1:Nbins;
    cosine             = ((cos(pi*(q)/(2*Nbins))).^cos_exp)';
    box                = (heavisideRF(q))';
    apodize_function   = cosine.*box;
    apo_NR             = signal_NR.*apodize_function;
    apo_RE             = signal_RE.*apodize_function;
    % FFT and phase
    NR_FT                               = fft(apo_NR,Nbins*zeropad_factor);
    RE_FT                               = fft(apo_RE,Nbins*zeropad_factor);
    if PowerSpectrum == 1
        dataStruct.PROC_2D_DATA{CurrDelay,1}   = abs((NR_FT + RE_FT).*exp(-1i*phases(CurrDelay)));
    else
        dataStruct.PROC_2D_DATA{CurrDelay,1}   = real((NR_FT + RE_FT).*exp(-1i*phases(CurrDelay)));
    end
    fittedPhase{CurrDelay,1}            = ones(Nbins)*phases(CurrDelay);
    
    % Pump correction and store in cell array
    % Write to dataStruct 
    dataStruct.PumpAxis{CurrDelay,1}       = PumpAxis';
    dataStruct.apodize_function    = apodize_function;
    dataStruct.t1delays            = t1delays; %
    dataStruct.binzero{CurrDelay,1}= 0; %
    dataStruct.binspecmax(CurrDelay,1) = dataStruct.binspecmax_pinhole;
%   dataStruct.FFT_ZPsig           = FFT_ZPsig;
%   dataStruct.phased_FFTZPsig     = phased_FFTZPsig;
% 	dataStruct.phased_FFTZPint		= phased_FFTZPint;
    dataStruct.fittedPhase			= fittedPhase;
%     dataStruct.phasepoints         = points;
    dataStruct.ZP_phase{CurrDelay,1}       = ones(Nbins*zeropad_factor).*phases(CurrDelay);
    dataStruct.phased_FFTZPint{CurrDelay,1}= dataStruct.powerspec;
    dataStruct.phase_coeff{CurrDelay,1}    = phases(CurrDelay);
%     dataStruct.apo_interferogram   = apo_interferogram;
%     dataStruct.apo_signal          = apo_signal;
    dataStruct.signal{CurrDelay,1} = signal_RE;
    dataStruct.SpecDiff            = 0;
end
% % contourf(PumpAxis,ProbeAxis,Power_spec')
% contourf(PumpAxis,ProbeAxis,PROC_2D_DATA{CurrDelay,1}')
% xlim([min(ProbeAxis) max(ProbeAxis)])
% ylim([min(ProbeAxis) max(ProbeAxis)])
% diagline = refline(1,0); diagline.Color = [0 0 0];
% colorbar

%% WRITE to dataStruct   
    dataStruct.isSimulation     = 0;
    dataStruct.Nspectra            = 1;
    dataStruct.Ndummies            = 1;
    dataStruct.Ndatastates         = 1;
    dataStruct.Nscans              = 1;