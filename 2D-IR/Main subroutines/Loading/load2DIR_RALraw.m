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
% Ricardo Fernandez-Teran / 09.11.2022 / v2.0a

%% HARDCODED Settings
DetSz       = [128 128];

%%% RAW data files (averaged)
SaveAvg     = 1;    % will save an averaged signal file
fixSign     = +1;   % In case I want to flip the signal again

%%% RAW data files (per cycle)
binShift    = -2;
signalSign  = -1; % the signal is upside down
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

% Try to get the number of runs with averaged files
% Ideally just 1 run but if name was not changed then there will be >1
runcyfil    = contains(datadir_fn(:),'run','IgnoreCase',1) & contains(datadir_fn(:),'avg','IgnoreCase',1) & contains(datadir_fn(:),'.csv','IgnoreCase',1);
theruncyfil = find(runcyfil);
RunCtr      = zeros(length(theruncyfil),1);
for i=1:length(theruncyfil)
    RunCtr(i) = sscanf(datadir_fn{theruncyfil(i)},'%*4c%d%*s.csv',1);
end
    
Nruns   = numel(unique(RunCtr));
thisRun = min(RunCtr);

if Nruns > 1
    warning('%i Runs Detected --- Will load run "%i". Move files to another folder to read other runs.',Nruns,thisRun);
end

% Check if there is a "home made" averaged file for this run, if there is, don't check again
AvgFilename  = [datadir filesep 'Run ' num2str(thisRun) '_AvgDiff.csv'];
AvgFileExist = exist(AvgFilename,'file');

% If there is NOT a "home made" averaged vile
if AvgFileExist ~= 2
    % Identify No. of Runs and No. of cycles
    % Ideally just 1 run but if name was not changed then there will be >1
    runcyfil    = contains(datadir_fn(:),'run','IgnoreCase',1) & contains(datadir_fn(:),'cycle','IgnoreCase',1) & contains(datadir_fn(:),'.csv','IgnoreCase',1);
    theruncyfil = find(runcyfil);
    RunCyCtr    = zeros(length(theruncyfil),2);
    for i=1:length(theruncyfil)
        RunCyCtr(i,:) = sscanf(datadir_fn{theruncyfil(i)},'%*4c%d%*6c%d.csv',[1 2]);
    end
        
    Nruns   = numel(unique(RunCyCtr(:,1)));
    thisRun = min(RunCyCtr(:,1));

    AvgFilename  = [datadir filesep 'Run ' num2str(thisRun) '_AvgDiff.csv'];
    if Nruns > 1
        warning('%i Runs Detected --- Will load run "%i". Move files to another folder to read other runs.',Nruns,thisRun);
    end
    
    % Figure out how many scans we will average. If Nscans = 0, something is wrong!
    Nscans  = max(RunCyCtr(RunCyCtr(:,1)==thisRun,2));
    
    % If something is missing, abort!
    if Nruns==0 || Nscans == 0
        error('Not a valid 2D-IR dataset: empty folder or corrupt data!')
    end
end

%% Check rotating frame and dt1
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
if AvgFileExist == 0 % Calculate from single scan data
    % Read the 1st scan
    % Create progress bar and clear error string
    if ShowWaitBar
        % dataStruct.WaitBar             = waitbar(0,'Loading data...');
        dataStruct.WBfigure                 = uifigure;
        dataStruct.WBfigure.Position(3:4)   = [405 175];
        dataStruct.WaitBar                  = uiprogressdlg(dataStruct.WBfigure,'Title','Loading 2D-IR data...','Message',['Loading RAW 2D-IR data... (Cycle 1 of ' num2str(Nscans) ')'],'Icon','info','ShowPercentage','on','Cancelable','on');
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
                dataStruct.WaitBar.Message  = ['Loading RAW 2D-IR data... (Cycle ' num2str(iScan) ' of ' num2str(Nscans) ')'];
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
                sig(:,:,p)     = real(log10(squeeze(signal4D(idxPrb,delID(m),:,p))).*phase_cy(p));
            end
            tmp = squeeze(sum(sig,3)).*signalSign.*1000;
            signal{m,k}  = tmp'; % CHECK WHETHER I NEED TO REMOVE LAST ELEMENT OR NOT?
            dummy_Onescell{m,k} = ones(Nbins,1);
            t1delays{m,k}       = [bins bins.*dt1];
        end
        est_probe{k} = (1:DetSz(k))';
    end

elseif AvgFileExist == 2
    %% If there is an averaged file, life is a bit easier...
    if ShowWaitBar
        % dataStruct.WaitBar             = waitbar(0,'Loading data...');
        dataStruct.WBfigure                 = uifigure;
        dataStruct.WBfigure.Position(3:4)   = [405 175];
        dataStruct.WaitBar                  = uiprogressdlg(dataStruct.WBfigure,'Title','Loading 2D-IR data...','Message','Loading RAW 2D-IR data... (Averaged data file)','Icon','info','ShowPercentage','on','Cancelable','on');
        dataStruct.WaitBar.Value            = 0;
        dataStruct.ErrorText.String         = "";
        drawnow;
        if dataStruct.WaitBar.CancelRequested
            delete(dataStruct.WBfigure);
            error('User aborted data loading!');
        end
    end
    
    % Read averaged file
    alldata = readmatrix(AvgFilename);
    
    if ShowWaitBar
        dataStruct.WaitBar.Value            = 0.5;
    end

    line1   = alldata(1,2:end);
    rawsig  = alldata(2:end,2:end);
    
    Npixels = size(rawsig,1);

    t2idx   = [0 find(diff(line1)<0)]+1;
    t2delays= line1(t2idx);
    t2delays= t2delays(:);
    Ndelays = length(t2delays);
    Nbins   = unique(diff(t2idx));
    
    % Generate bins and delays axes
    bins        = (1:Nbins)'-1;
    t1delays    = cell(Ndelays,2);

    % Do signal calc.
    signal      = cell(Ndelays,1);
    est_probe   = cell(2,1);
    cmprobe     = cell(1,2);
    dummy_Onescell = cell(Ndelays,2);

    Ndet = length(DetSz);
    for k=1:Ndet
        idxPrb = (1:DetSz(k)) + sum(DetSz(1:k-1));
        for m=1:Ndelays
            idxBins             = (m-1)*Nbins + (1:Nbins);
            signal{m,k}         = rawsig(idxPrb,idxBins)'.*fixSign; % If I want to flip the signal again
            dummy_Onescell{m,k} = ones(Nbins,1);
            t1delays{m,k}       = [bins bins.*dt1];
        end
        est_probe{k} = (1:DetSz(k))';
    end
end

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
    dataStruct.signal        = signal; % cell array of [Nbins x Npix(per det)]
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


if AvgFileExist == 0 && SaveAvg == 1
    avgDat = [];
    for m=1:Ndelays
        avgDat = [avgDat [[t2delays(m), t2delays(m)+10e6.*((1:Nbins-1)-1)]; [signal{m,1} signal{m,2}]']]; %#ok<*AGROW> 
    end
    avgDat = [(0:2*Npixels)' avgDat];
    tic
        writematrix(avgDat,AvgFilename);
    tw=toc;
    fprintf('Wrote averaged 2DIR data file in %f s\n\tFilename: \t%s\n\n',tw,AvgFilename);
end