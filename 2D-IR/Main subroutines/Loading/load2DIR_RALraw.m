function  dataStruct = load2DIR_RALraw(dataStruct,varargin)

% Description: This function loads RAW 2DIR data from the ULTRA LIFEtime experiment at the
% Central Laser Facility of the Rutherford Appleton Laboratory.
%
% Usage: dataStruct = load2DIR_RALraw(dataStruct)
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
% Ricardo Fernandez-Teran / 05.11.2022 / v1.6b

%% DEBUG
% rootdir         = 'D:\GoogleDrive\RESULTS\From James\2DIR DATA';
% datafilename    = 'RAL LIFEtime - LT5 PTZCCC13C13NAP  2DIR 22 fs';
% ShowWaitBar     = true;
% 

%% HARDCODED Settings
binShift    = -2;
DetSz       = [128 128];
signalSign  = -1; % the signal is upside dowm

Nphases     = 4;
phase_cy    = [+1 -1 +1 -1]; % SIGNAL
% phase_cy    = [+1 +1 +1 +1]; % SCATTER
% phase_cy    = [+1 +1 -1 -1]; % SCATTER?


%% READ from dataStruct
if isempty(varargin) || varargin{1} == 1
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
fn          = fn(contains(fn,'Avg Diff','ignorecase',1));

% Identify No. of Runs (!)
% Ideally just 1 run but if name was not changed then there will be >1
Nruns = sum(contains(datadir_fn(:),'run','IgnoreCase',1) & contains(datadir_fn(:),'avg','IgnoreCase',1));
runID = [];
for i = find(contains(datadir_fn(:),'avg','IgnoreCase',1))'
    [~,fnRuns,~]    =  fileparts(datadir_fn{i});
    runID = [runID sscanf(fnRuns,'Run %i %s')];
end
thisRun = runID(1);

if Nruns > 1
    warning('%i Runs Detected --- Will load run "%i". Move files to another folder to read other runs.',Nruns,thisRun);
end

% Figure out how many scans we will average. If Nscans = 0, something is wrong!
Nscans      = sum(contains(datadir_fn(:),['run ' num2str(thisRun)],'IgnoreCase',1) & contains(datadir_fn(:),'cycle','IgnoreCase',1) & contains(datadir_fn(:),'.csv','IgnoreCase',1));

if ~isempty(fn) || isempty(fnRuns) || Nscans == 0
% If there is an w0 file in the current ROOTDIR, use it to set the rotating frame frequency.
% Otherwise, plot relative to w0
if exist([datadir filesep 'w0.csv'],'file') == 2
	rotframe = readmatrix([datadir filesep 'w0.csv']);
    w0  = rotframe(1);
    dt1 = rotframe(2);
elseif varargin{2} == 0
    w0  = 0;
    dt1 = 1;
else
    prompt = {'Enter rotating frame frequency (\omega_{r}, in cm^{-1}):','Enter t_{1} time step (in fs):'};
    dlgtitle = 'Rotating Frame Settings';
    definput = {'1200','15'};
    dims = [1 40];
    opts.Interpreter = 'tex';
    answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
    
    if isempty(answer)
        w0  = 0;
        dt1 = 1;
    else
        w0  = str2double(answer{1});
        dt1 = str2double(answer{2});
        try
            writematrix([w0;dt1],[datadir filesep 'w0.csv']);
        catch
            warning('Error writing rotating frame info file to default path');
            fp = uigetdir('Where to save rotating frame file?','w0.csv');
            if fp == 0
                warning('Nowhere to save the rotating frame file :(');
            else
                writematrix([w0;dt1],[fp filesep 'w0.csv']);
            end          
%             warndlg('Error writing rotating frame info file','Error writing file');
        end
    end
end

%% Read data
% Read the 1st scan
% Create progress bar and clear error string
if ShowWaitBar
    % dataStruct.WaitBar             = waitbar(0,'Loading data...');
    dataStruct.WBfigure                 = uifigure;
    dataStruct.WBfigure.Position(3:4)   = [405 175];
    dataStruct.WaitBar                  = uiprogressdlg(dataStruct.WBfigure,'Title','Loading 2D-IR data...','Message',['Loading RAW 2D-IR data... (Scan 1 of ' num2str(Nscans) ')'],'Icon','info','ShowPercentage','on','Cancelable','on');
    dataStruct.WaitBar.Value            = 0.5./Nscans;
    dataStruct.ErrorText.String         = "";
    drawnow;
    if dataStruct.WaitBar.CancelRequested
        delete(dataStruct.WBfigure);
        error('User aborted data loading!');
    end
end

cycleData                   = readmatrix([datadir filesep 'Run ' num2str(thisRun) ' Cycle1.csv']);
if ShowWaitBar
    dataStruct.WaitBar.Value    = (0.5*1./Nscans);
    drawnow;
end

% Read the rest of the scans and average them
if Nscans > 1
    for iScan=2:Nscans
        if ShowWaitBar
            % Update the Wait Bar
            dataStruct.WaitBar.Value    = (0.5*iScan./Nscans);
            dataStruct.WaitBar.Message  = ['Loading RAW 2D-IR data... (Scan ' num2str(iScan) ' of ' num2str(Nscans) ')'];
            if dataStruct.WaitBar.CancelRequested
                delete(dataStruct.WBfigure);
                error('User aborted data loading!');
            end
            drawnow;
        end
        cycleData   = cycleData + readmatrix([datadir filesep 'Run ' num2str(thisRun) ' Cycle' num2str(iScan) '.csv']);
    end
end
% Divide by Nscans to calculate average
cycleData = cycleData ./ Nscans;

%% Figure out how many masks and t2 delays we have

cyPixData   = cycleData(2:end,:);
line1       = cycleData(1,:)';
line1(isnan(line1)) = 0;

Nwfm        = find(diff(line1) < 0,1)-1;
Ndelays     = (length(line1)-1)./(2.*Nwfm-1);
Nbins       = Nwfm./Nphases;

delidx      = 2 + [0 (1:Ndelays-1).*(2.*Nwfm-1)]' + 1; %#ok<*BDSCI>             % adapted from ULTRAvew for 2DIR (v1.8)

t1ranges    = [delidx delidx+Nwfm-1]-1;

t2delays    = line1(delidx);
%%

if mod(Nbins,1) ~= 0 % Check for non-integer number of bins
    error('Inconsistent number of bins for phase cycling, please check!');
end

Ndelays     = length(t2delays);
Npixels     = size(cyPixData,1);

%% Signal calc
signal4D    = zeros(Npixels,Ndelays,Nbins,Nphases);

for i=1:Ndelays
    % index of all masks for this delay
    maskIndex       = circshift(t1ranges(i,1):t1ranges(i,2),binShift); % Shift the mask sequence here (!)
    modPhase        = [0 1:Nphases-1];
    for j=1:Nphases
        idx                 = maskIndex(mod(maskIndex - t1ranges(i,1),Nphases)==modPhase(j));
        signal4D(:,i,:,j)   = cyPixData(:,idx); 
    end
end

%% Calculate Actual Signal from Phase Cycling [ASSUMING 4 FRAME, NO CHOP: p1{0,0,pi,pi} p2{0,pi,0,pi}]
% Signal = + A - B + C - D

signal      = cell(Ndelays,1);
est_probe   = cell(2,1);
cmprobe     = cell(1,2);
dummy_Onescell = cell(Ndelays,2);

% Sort t2 delays
[t2delays,delID] = sort(t2delays);

% Generate bins and t1 axis
% Nbins       = size(tmp,2) - 0; % CHECK WHETHER I NEED TO REMOVE LAST ELEMENT OR NOT?
bins        = (1:Nbins)'-1;
t1delays    = cell(Ndelays,2);

% Do signal calc.
Ndet = length(DetSz);
for k=1:Ndet
    idxPrb = (1:DetSz(k)) + sum(DetSz(1:k-1));
    Npixels= length(idxPrb);
    sig    = zeros(Npixels,Nbins,Nphases);
    for m=1:Ndelays
        for p=1:Nphases
            sig(:,:,p)     = log10(squeeze(signal4D(idxPrb,delID(m),:,p))).*phase_cy(p);
        end
        tmp = squeeze(sum(sig,3)).*signalSign.*1000;
        signal{m,k}  = tmp'; % CHECK WHETHER I NEED TO REMOVE LAST ELEMENT OR NOT?
        dummy_Onescell{m,k} = ones(Nbins,1);
        t1delays{m,k}       = [bins bins.*dt1];
    end
    est_probe{k} = (1:DetSz(k))';
end

% plot(signal{10,2}(:,59:68)); xlim([-1,10])

%% WRITE to dataStruct (Load)
    dataStruct.Gratings      = 0;
    dataStruct.CWL           = 0;
    dataStruct.est_probe     = est_probe; % Estimated probe from spectrograph info
    dataStruct.DetSz         = DetSz;

    dataStruct.isSimulation  = 0;
    dataStruct.isShaper      = 1;
    dataStruct.cmprobe       = cmprobe;
    dataStruct.bins          = bins;
    dataStruct.t2delays      = t2delays; % Delays already in ps !!!
    dataStruct.Ndelays       = Ndelays;
    dataStruct.signal        = signal;
    dataStruct.Nspectra      = 2;
    dataStruct.Ndummies      = 1;
    dataStruct.Ndatastates   = 1;
    dataStruct.Nbins         = Nbins;
    dataStruct.Nslowmod      = 1;
    dataStruct.t1delays      = t1delays; % in fs!
    dataStruct.Nscans        = 1;
    dataStruct.datatype      = 'TimeFreq';
    dataStruct.interferogram = dummy_Onescell;
    dataStruct.w0            = w0;
    dataStruct.dt1           = dt1;
% Clear 2DGC fit results
dataStruct.FitResults    = [];
dataStruct.t2_startFit   = [];
dataStruct.FitInput      = [];
else
    error('Not a valid 2D-IR dataset: empty folder or corrupt data')
end