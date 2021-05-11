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
% Ricardo Fernandez-Teran / 11.05.2021 / v0.9a

%% DEBUG
% rootdir         = 'D:\GoogleDrive\RESULTS\From James';
% datafilename    = 'PTZCCC13C13NAP LT5';
% ShowWaitBar     = true;

%% HARDCODED Settings

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

[~,fn,~]    = fileparts({datadir_fn{contains(datadir_fn(:),'.csv','IgnoreCase',1)}});

Ndelays     = length(fn);
fn_cell     = cell(Ndelays,1);
t2delays_u    = zeros(Ndelays,1);

for i=1:Ndelays
    fn_cell{i}  = strsplit(fn{i},'_');
    t2delays_u(i) = str2double(erase(fn_cell{i}{end-1},'ps'));
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

rawdata   = readmatrix([datadir filesep fn{1} '.csv']);
    
Nbins     = size(rawdata,2)-1;
bins      = (1:Nbins)';

cmprobe   = rawdata(2:end,1);  % (?)

% [cmprobe,idProbe,~] = unique(cmprobe);
Npixels   = length(cmprobe);

lowpart   = 1:(Npixels/2);
highpart  = (Npixels/2:Npixels-1)+1;

ProbeCut  = 'other';

switch ProbeCut
    case 'low'
        startlow    = 1;
        endlow      = length(cmprobe(cmprobe(lowpart) < min(cmprobe(highpart))));
        starthigh   = Npixels/2+1;
        endhigh     = Npixels;
    case 'high'
        startlow    = 1;
        endlow      = Npixels/2;
        starthigh   = Npixels/2 + find(cmprobe(highpart) > max(cmprobe(lowpart)),1);
        endhigh     = Npixels;
    case 'other'
        startlow    = 1;
        endlow      = Npixels/2 - 20;
        starthigh   = Npixels/2 + find(cmprobe(highpart) > max(cmprobe(startlow:endlow)),1);
        endhigh     = Npixels;
end
idProbe = [startlow:endlow starthigh:endhigh];
Npixels = length(idProbe);
cmprobe = cmprobe(idProbe);
% idProbe   = 1:Npixels;

signal          = cell(Ndelays,1);
t1delays        = cell(Ndelays,1);
dummy_cell      = cell(Ndelays,1);
dummy_Onescell  = cell(Ndelays,1);
PumpAxis        = cell(Ndelays,1);
PROC_2D_DATA    = cell(Ndelays,1);

for i=1:Ndelays
    dummy_cell{i,1}     = 0;
    dummy_Onescell{i,1} = ones(Nbins,1);
    signal{i,1}         = zeros(Nbins,Npixels);
    t1delays{i,1}       = bins;

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
    
    rawdata = readmatrix([datadir filesep fn{idx(i)} '.csv']);
    PumpAxis{i,1}       = rawdata(1,2:end)';
    PROC_2D_DATA{i,1}   = 1000*rawdata(idProbe+1,2:end)';
end
%% WRITE to dataStruct (Load)
    dataStruct.isSimulation  = 0;
    dataStruct.isShaper      = 1;
    dataStruct.cmprobe       = cmprobe;
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
    dataStruct.datatype      = 'TimeFreq';
    dataStruct.interferogram = dummy_Onescell;

%% WRITE to dataStruct (Process 2D-IR)
    dataStruct.ProbeAxis           = cmprobe;
    dataStruct.freq_fit            = [];
    dataStruct.scattering_maxima   = [];
    dataStruct.PumpAxis            = PumpAxis;
    
    dataStruct.t1delays            = t1delays;
    
    dataStruct.binzero             = dummy_cell;
    dataStruct.binspecmax          = ones(Ndelays,1);
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