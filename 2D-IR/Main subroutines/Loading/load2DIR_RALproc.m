function  dataStruct = load2DIR_RALproc(dataStruct,varargin)

% Description: This function loads pre-processed 2DIR data from the ULTRA LIFEtime experiment at the
% Central Laser Facility of the Rutherford Appelton Laboratory.
%
% Usage: dataStruct = load2DIRsimu(dataStruct)
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
%     count             (Cell array)
%     probe             (Cell array)
%     reference         (Cell array)
%     interferogram     (Cell array)
%     signal			(Cell array)
%
% Ricardo Fernandez-Teran / 22.04.2022 / v1.1b

%% HARDCODED Settings
DetSz = [128 128]; % Probe 1 / Probe 2
dt1 = 1; % Since the data is already in frequency domain

SIGN = +1; % In case I want to invert the data for whatever reason (e.g. wrong signal calculation)

%% READ from dataStruct
if isempty(varargin)
    ShowWaitBar = true;
else
    ShowWaitBar = false;
end

datafilename    = dataStruct.datafilename;
rootdir         = dataStruct.rootdir;

%% Load all the necessary files after checking that they exist
datadir     = [rootdir filesep datafilename];
datadir_fl  = dir(datadir);
datadir_fn  = {datadir_fl.name};

[~,fn,~]    = fileparts({datadir_fn{contains(datadir_fn(:),'2DIR.csv','IgnoreCase',1)}});

Ndelays     = length(fn);
fn_cell     = cell(Ndelays,1);
t2delays_u  = zeros(Ndelays,1);

for m=1:Ndelays
    fn_cell{m}    = strsplit(fn{m},'_');
    t2delays_u(m) = str2double(erase(fn_cell{m}{end-1},'ps'));
end

[t2delays,idx] = sort(t2delays_u); % Delays are already in ps

% Create progress bar and clear error string
if ShowWaitBar == 1
    progress = 0;
    % dataStruct.WaitBar             = waitbar(0,'Loading data...');
    dataStruct.WBfigure                 = uifigure;
    dataStruct.WBfigure.Position(3:4)   = [405 175];
    dataStruct.WaitBar                  = uiprogressdlg(dataStruct.WBfigure,'Title','2D-IR data processing...','Message','Loading data...','Icon','info','ShowPercentage','on','Cancelable','on');
    dataStruct.ErrorText.String         = "";
end


if Ndelays > 0
rawdata = readmatrix([datadir filesep fn{1} '.csv']);
    
Nbins   = size(rawdata,2)-1;
bins    = (1:Nbins)';
w3      = rawdata(2:end,1);

ProbeAxis       = cell(2,1);
signal          = cell(Ndelays,2);
t1delays        = cell(Ndelays,2);
dummy_cell      = cell(Ndelays,1);
dummy_Onescell  = cell(Ndelays,1);
PumpAxis        = cell(Ndelays,1);
PROC_2D_DATA    = cell(Ndelays,1);

for m=1:Ndelays 
    if ShowWaitBar == 1
        % Update the Wait Bar
        progress = progress + 1;
        dataStruct.WaitBar.Value    = (progress/(Ndelays));
        dataStruct.WaitBar.Message  = ['Processing data... (' num2str(progress) ' of ' num2str(Ndelays) ')'];
        if dataStruct.WaitBar.CancelRequested
            delete(dataStruct.WBfigure);
            error('User aborted loading the data!');
        end
        drawnow;
    end
    
    rawdata = readmatrix([datadir filesep fn{idx(m)} '.csv']);
    for k=1:2 % Ndetectors
        ProbeIdx = (1:DetSz(k)) + sum(DetSz(1:k-1));
        dummy_cell{m,k}     = 0;
        dummy_Onescell{m,k} = ones(Nbins,1);
        binzero{m,k}        = 1; % Because we are using the shaper
        binspecmax(m,k)     = 1; %#ok<*AGROW> 
        
        signal{m,k}         = zeros(Nbins,DetSz(k));
        t1delays{m,k}       = [bins bins.*dt1];
        PumpAxis{m,k}       = rawdata(1,2:end)';
        PROC_2D_DATA{m,k}   = SIGN*1000*rawdata(ProbeIdx,2:end)';
    end
end
for k=1:2
    ProbeIdx = (1:DetSz(k)) + sum(DetSz(1:k-1));
    ProbeAxis{k} = w3(ProbeIdx);
end
%% WRITE to dataStruct (Load)
    dataStruct.isSimulation  = 0;
    dataStruct.isShaper      = 1;
    dataStruct.cmprobe       = ProbeAxis;
    dataStruct.bins          = bins;
    dataStruct.t2delays      = t2delays; % Delays already in ps !!!
    dataStruct.Ndelays       = Ndelays;
    dataStruct.signal        = signal;
    dataStruct.Nspectra      = 1;
    dataStruct.Ndummies      = 1;
    dataStruct.Ndatastates   = 1;
    dataStruct.Nbins         = Nbins;
    dataStruct.Nslowmod      = 1;
    dataStruct.t1delays      = t1delays; % in fs!
    dataStruct.Nscans        = 1;
    dataStruct.datatype      = 'FreqFreq';
    dataStruct.interferogram = dummy_Onescell;

%% WRITE to dataStruct (Process 2D-IR)
    dataStruct.ProbeAxis           = ProbeAxis;
    dataStruct.freq_fit            = dummy_Onescell;
    dataStruct.scattering_maxima   = dummy_Onescell;
    dataStruct.PumpAxis            = PumpAxis;
    
    dataStruct.t1delays            = t1delays;
    
    dataStruct.binzero             = binzero;
    dataStruct.binspecmax          = binspecmax;
    dataStruct.apodize_function    = dummy_Onescell;
    dataStruct.FFT_ZPsig           = [];
    dataStruct.phased_FFTZPsig     = dummy_Onescell;
	dataStruct.phased_FFTZPint     = dummy_Onescell;
    dataStruct.fittedPhase         = dummy_Onescell;
    dataStruct.phasepoints         = dummy_Onescell;
    dataStruct.ZP_phase            = dummy_cell;
    dataStruct.phase_coeff         = dummy_cell;
    dataStruct.apo_interferogram   = [];
    dataStruct.apo_signal          = signal;
    dataStruct.PROC_2D_DATA        = PROC_2D_DATA;
    dataStruct.SpecDiff            = 0;
    
% Clear 2DGC fit results
dataStruct.FitResults    = [];
dataStruct.t2_startFit   = [];
dataStruct.FitInput      = [];

else
    error('Not a valid 2D-IR dataset: empty folder or corrupt data')
end