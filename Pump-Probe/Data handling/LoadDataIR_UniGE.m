function dataStruct = LoadDataIR_UniGE(dataStruct,recalcScans)
%% READ from dataStruct
rootdir         = dataStruct.rootdir;
datafilename    = dataStruct.datafilename;

% Decide whether to flip the sign of the signal
sign=+1;

%% Check that all files exist and read them
% Get files
datadir     = [rootdir filesep datafilename];
datadir_fl  = dir(datadir);
datadir_fn  = {datadir_fl.name};

scanfiles   = datadir_fn(contains(datadir_fn(:),'.top','IgnoreCase',true));

Nscans      = length(scanfiles);

[~,fn,~]    = fileparts(scanfiles);

fullName    = [rootdir filesep datafilename filesep fn{1}];
% If there is a calibrated WL file in the current EXPDIR, use it
wavenumberfile  = 'CalibratedProbe.csv';
if exist([datadir filesep wavenumberfile],'file') == 2 % it's inside the data folder
    probe_calib = 2;
    tmp_probe = readmatrix([datadir filesep wavenumberfile]);
elseif exist([rootdir filesep wavenumberfile],'file') == 2 % it's outside
    probe_calib = 2;
    tmp_probe   = readmatrix([rootdir filesep wavenumberfile]);
else % if it doesn't exist, takewhat we get from the spectrograph
    probe_calib = 0;
    tmp_probe   = readmatrix([fullName '.wvl'],'FileType','text','Delimiter','\t');
    tmp_probe   = 1e7./tmp_probe;
end
cmprobe = {tmp_probe(:)};

% Read delays
delays = readmatrix([fullName '.stg1'],'FileType','text','Delimiter','\t');

% Read first scan data and initialise variables
scandata{1}(:,:,1)  = -1e3.*sign.*readmatrix([fullName '.top'],'FileType','text','Delimiter','\t');

Ndelays = size(scandata{1},1);
if Ndelays ~= length(delays)
    error('Uh oh! Inconsistent number of delays. Check the data and try again!')
end

Npixels = size(scandata{1},2);
if Npixels ~= length(cmprobe{1})
    error('Uh oh! Inconsistent number of pixels. Check the data and try again!')
end

rawsignal{1} = zeros(Ndelays,Npixels);
noise{1}     = zeros(Ndelays,Npixels);

%% Read single scan data 
if Nscans >= 1
    scandata{1}         = zeros([size(rawsignal{1}) Nscans]);
    scan_err{1}         = zeros([size(rawsignal{1}) Nscans]);
    for s=2:Nscans
        fullName    = [rootdir filesep datafilename filesep fn{s}];
        scandata{1}(:,:,s)  = -1e3.*sign.*readmatrix([fullName '.top'],'FileType','text','Delimiter','\t');
    end
    noise{1}            = std(scandata{1},0,3); % in mOD
    rawsignal{1}        = mean(scandata{1},3); % in mOD
else
    Nscans      = NaN;
    rawsignal{1}= scandata{1}(:,:,1);
    noise{1}    = zeros(size(rawsignal{1}));
    scandata    = [];
    scan_err    = [];
end

if recalcScans == 1
    rawsignal       = dataStruct.rawsignal;
    noise           = dataStruct.noise;
    cmprobe         = dataStruct.cmprobe;
    delays          = dataStruct.delays;
    Nscans          = 0;
end

%% Read the plot ranges
mintime     = min(delays);
maxtime     = max(delays);
minabs      = min(rawsignal{1}(:));
maxabs      = max(rawsignal{1}(:));
zminmax     = round(max([abs(minabs) abs(maxabs)]),3);
minwl       = min(cmprobe{1});
maxwl       = max(cmprobe{1});
Ncontours   = 40; % 40 contours by default is OK
plotranges{1} = [mintime maxtime minwl maxwl minabs maxabs Ncontours];

%% WRITE to dataStruct
% Write the main variables
dataStruct.delays      = delays;
dataStruct.cmprobe     = cmprobe;
dataStruct.rawsignal   = rawsignal;
dataStruct.noise       = noise;
dataStruct.plotranges  = plotranges;
dataStruct.scandata    = scandata;
dataStruct.scanNoise   = scan_err;

% Calculate noise statistics
dataStruct.AvgNoise    = mean(noise{1}(:),'omitnan');
dataStruct.MaxNoise    = max(noise{1}(:));
dataStruct.SNR         = abs(round(zminmax/dataStruct.AvgNoise,3));

% Number of scans
dataStruct.Nscans      = Nscans;

% Background Subtraction
if dataStruct.recalcBkg == 0
    % Subtract the first negative delay from all the dataset by default - otherwise take the inputs
    dataStruct.mintimeBkg = dataStruct.delays(1);
    dataStruct.maxtimeBkg = dataStruct.delays(3);
end
    Idx = findClosestId2Val(dataStruct.delays,[dataStruct.mintimeBkg dataStruct.maxtimeBkg]);
    % Do the background subtraction and change status in handles.rawcorr
    if Idx(2)==1
        dataStruct.bkg{1} = dataStruct.rawsignal{1}(1,:);
    else
        dataStruct.bkg{1} = mean(dataStruct.rawsignal{1}(Idx(1):Idx(2),:),'omitnan');
    end
dataStruct.corrdata{1}  = dataStruct.rawsignal{1} - dataStruct.bkg{1};

dataStruct.timescale = 'ps';
fclose('all');