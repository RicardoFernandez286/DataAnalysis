function [dataStruct,isTOPAS,isDelay,TOPASwn] = LoadDataIR_ULTRA(dataStruct,recalcScans)
%% READ from dataStruct
% Get stuff from dataStruct 
rootdir     = dataStruct.rootdir;
datafilename= dataStruct.datafilename;

NActiveDet = 2;

DetSz   = [128 128]; % Sizes of [probe1 probe2 ref1 ref2] in pixels
SortScans   = 1;

%% Get files
datadir     = [rootdir filesep datafilename];
datadir_fl  = dir(datadir);
datadir_fn  = {datadir_fl.name};

% Identify No. of Runs (!)
% Ideally just 1 run but if name was not changed then there will be >1
Nruns = sum(contains(datadir_fn(:),'run','IgnoreCase',1) & contains(datadir_fn(:),'avg diff','IgnoreCase',1));
runID = [];
for i = find(contains(datadir_fn(:),'avg diff','IgnoreCase',1))'
    [~,fn1,~]    =  fileparts(datadir_fn{i});
    runID = [runID sscanf(fn1,'Run %i %s')];
end
thisRun = runID(1);

if Nruns > 1
    warning('%i Runs Detected --- Will load run "%i". Move files to another folder to read other runs.',Nruns,thisRun);
end

Nscans      = sum(contains(datadir_fn(:),['run ' num2str(thisRun)],'IgnoreCase',1) & contains(datadir_fn(:),'cycle','IgnoreCase',1) & contains(datadir_fn(:),'.csv','IgnoreCase',1));

% [~,fn,~]    = fileparts(datadir_fn{contains(datadir_fn(:),'run 0','IgnoreCase',1) & contains(datadir_fn(:),'avg','IgnoreCase',1)});
% 


%% Probe Calibration
% If there is a calibrated WL file in the current EXPDIR, use it
wavenumberfile  = 'CalibratedProbe.csv';

if exist([datadir filesep wavenumberfile],'file') == 2
    probe_calib = 2;
    tmp_probe = readmatrix([datadir filesep wavenumberfile]);
elseif exist([rootdir filesep wavenumberfile],'file') == 2
    probe_calib = 2;
    tmp_probe   = readmatrix([rootdir filesep wavenumberfile]);
else  
    probe_calib = 0;
    est_probe   = 1:256;
end

switch probe_calib
    case 0
        cmprobe         = est_probe;
    case 2
        cmprobe{1}      = tmp_probe(1:128);
        if NActiveDet == 2
            cmprobe{2}  = tmp_probe(129:end);
        end
end

%% Read the single scan files
% Read first scan to get array sizes
iScan = 1;
    [~,fn,~]    = fileparts(datadir_fn{contains(datadir_fn(:),['run ' num2str(thisRun)],'IgnoreCase',1) & contains(datadir_fn(:),['Cycle' num2str(iScan) '.csv'],'IgnoreCase',1)});
    filename    = [rootdir filesep datafilename filesep fn];
    scanData    = readmatrix([filename '.csv'],'FileType','delimited');
    fileID = fopen([filename '.csv'],'r');
    L1 = fgetl(fileID); L1 = strsplit(L1,','); L1=string(L1{1});
    L2 = fgetl(fileID); L2 = strsplit(L2,','); L2=string(L2{1});
    L3 = fgetl(fileID); L3 = strsplit(L3,','); L3=string(L3{1});
    fclose(fileID);
    
    isTOPAS = contains([L1;L2;L3],'TOPAS','IgnoreCase',true);
    isDelay = contains([L1;L2;L3],'Delay','IgnoreCase',true);
    
    % If there is a TOPAS scan inside
    if sum(isTOPAS) > 0
        TOPASwn = scanData(isTOPAS,2:end)';
        TOPASwn = unique(TOPASwn(~isnan(TOPASwn)),'rows');
        N_TOPAS = numel(TOPASwn);
    else
        N_TOPAS = 1;
        TOPASwn = NaN;
    end
    
    % If there are delay scans inside
    if sum(isDelay) > 0
        delays  = scanData(isDelay,2:end)';
        delays  = unique(delays(~isnan(delays)),'rows');
        Ndelays = numel(delays);
    else
        Ndelays = 1;
    end

    % Decide how the data is going to be sorted. It contains
    % UV-IR/IR-IR/t2DIR [3 spectra]
    % These are in different files, so I would sort it as [UV|IR|t2D] for
    % each TOPAS position and detector. 
    A=0;
    %%% CONTINUE HERE

    % Every third column contains a valid delay
    scan1data   = scanData(2:end,4:3:end)';
    scan1data(abs(scan1data)==Inf) = NaN;

    
    rawdata3D   = zeros([size(scan1data) Nscans]);
    rawdata3D(:,:,1)  = scan1data;

if Nscans > 1
    for iScan=2:Nscans
        [~,fn,~]    = fileparts(datadir_fn{contains(datadir_fn(:),['run ' num2str(thisRun)],'IgnoreCase',1) & contains(datadir_fn(:),['Cycle' num2str(iScan) '.csv'],'IgnoreCase',1)});
        filename    = [rootdir filesep datafilename filesep fn];
        scanData    = readmatrix([filename '.csv'],'FileType','delimited');       
        scanData(abs(scanData)==Inf) = NaN;
        rawdata3D(:,:,iScan) = scanData(2:end,4:3:end)';
    end
end

rawsignal       = cell(2,1);
plotranges      = cell(2,1);
noise           = cell(2,1);
scandata        = cell(2,1);
scannoise       = cell(2,1);

if recalcScans == 1
    rawsignal   = dataStruct.rawsignal;
    noise       = dataStruct.noise;
    cmprobe     = dataStruct.cmprobe;
    delays      = dataStruct.delays;
    Nscans      = 0;
else
    %% Split Data into two detectors
    [delays,delID]  = sort(delays);

    for j=1:2
        ProbeIdx = (1:DetSz(j)) + sum(DetSz(1:j-1));
        if SortScans==1
            scandata{j}  = rawdata3D(delID,ProbeIdx,:).*1000; % to mOD
            rawsignal{j} = mean(scandata{j},3,'omitnan'); % mean across scans
            scannoise{j} = zeros(size(rawsignal{j})); % no single scan noise information
            noise{j}     = std(scandata{j},0,3,'omitnan').*1000; % to mOD
        else
            Nscans      = NaN;
            scandata    = [];
            scannoise   = [];
            noise{j}    = zeros(size(rawsignal{j})); % no noise information
        end
    end
   
end

% Set timescale
if max(delays) > 5000
    delays = delays./1000;
    dataStruct.timescale    = 'ns';
else
    dataStruct.timescale    = 'ps';
end

for j=1:NActiveDet
    % Noise and ranges
    
    % Read the plot ranges
    mintime         = min(delays);
    maxtime         = max(delays);
    minabs          = min(rawsignal{j},[],'all','omitnan');
    maxabs          = max(rawsignal{j},[],'all','omitnan');
    zminmax         = round(max([abs(minabs) abs(maxabs)]),3);
    minwl           = min(cmprobe{j});
    maxwl           = max(cmprobe{j});
    Ncontours       = 40; % 40 contours by default is OK
    plotranges{j}   = [mintime maxtime minwl maxwl minabs maxabs Ncontours];
end

%% WRITE to dataStruct
% Write the main variables
dataStruct.delays       = delays;
dataStruct.cmprobe      = cmprobe;
dataStruct.rawsignal    = rawsignal;
dataStruct.noise        = noise;
dataStruct.plotranges   = plotranges;

% Calculate noise statistics
dataStruct.AvgNoise     = mean(noise{2}(:));
dataStruct.MaxNoise     = max(noise{1}(:),[],'all','omitnan');
dataStruct.SNR          = abs(round(zminmax/dataStruct.AvgNoise,3));

% Number of scans
dataStruct.Nscans       = Nscans;
dataStruct.scandata     = scandata;
dataStruct.scanNoise    = scandata;

%% Background correction
if recalcScans == 1
    dataStruct.recalcBkg = 0;
end

for j=1:NActiveDet
    % Background Subtraction
    if dataStruct.recalcBkg == 0
        % Subtract the first negative delay from all the dataset by default - otherwise take the inputs
        dataStruct.mintimeBkg = dataStruct.delays(1);
        dataStruct.maxtimeBkg = dataStruct.delays(1);
    end
        Idx = findClosestId2Val(dataStruct.delays,[dataStruct.mintimeBkg dataStruct.maxtimeBkg]);
        % Do the background subtraction and change status in handles.rawcorr
        if Idx(2)==1
            dataStruct.bkg{j} = dataStruct.rawsignal{j}(1,:);
        else
            dataStruct.bkg{j} = mean(dataStruct.rawsignal{j}(Idx(1):Idx(2),:));
        end
    dataStruct.corrdata{j} = dataStruct.rawsignal{j} - dataStruct.bkg{j};
end

