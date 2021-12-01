function dataStruct = LoadDataIR_UoS(dataStruct,varargin)
%% READ from dataStruct
% Get stuff from dataStruct
rootdir         = dataStruct.rootdir;
datafilename    = dataStruct.datafilename;

DetSz           = [96 96 32 32]; % Sizes of [probe1 probe2 ref1 ref2] in pixels

% Get files
datadir     = [rootdir filesep datafilename];
datadir_fl  = dir(datadir);
datadir_fn  = {datadir_fl.name};

[~,fn,~]    = fileparts(datadir_fn{contains(datadir_fn(:),'.LG','IgnoreCase',1)});

filename    = [rootdir filesep datafilename filesep fn];

% Read Spectrograph Information
text        = readlines([filename '.LG']);
g1_st       = strsplit(text{75},' '); % Det186 is 1st
w1_st       = strsplit(text{77},' ');
g2_st       = strsplit(text{67},' '); % Det185 is 2nd
w2_st       = strsplit(text{69},' ');

Gratings(1) = str2double(g1_st{end});
CWL(1)      = str2double(w1_st{end});
Gratings(2) = str2double(g2_st{end});
CWL(2)      = str2double(w2_st{end});

probe_calib = 0;
est_probe   = cell(1,2);

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

% Estimate probe axis, arbitrary calibration factors
for i=1:2
    switch Gratings(i)
        case 0 % 120 l/mm
            ppnm    = 0.168;
        case 1 % 100 l/mm
            ppnm    = 0.14;
        case 2 % 50 l/mm
            ppnm    = 0.068;
    end
    pix = flip(1:DetSz(i));
    est_probe{i} = 1e7./(pix./ppnm + CWL(i)+60 - DetSz(i)./2./ppnm);
end

%% Read the file
rawdata         = readmatrix([filename '.2D'],'FileType','delimited');
delays          = readmatrix([filename '.DT'],'FileType','delimited'); % Delays originally in microseconds, convert to seconds

switch probe_calib
    case 0
        cmprobe         = est_probe;
    case 2
        cmprobe{1}      = tmp_probe(1:96);
        cmprobe{2}      = tmp_probe(97:end);
end

Ndelays         = length(delays);
rawsignal       = cell(2,1);
%% Split Data into two detectors
for j=1:2
    idx = (1:DetSz(j)) + sum(DetSz(1:j-1));
    rawsignal{j} = rawdata(:,idx);

    % Noise and ranges
    noise           = zeros(size(rawsignal));
    % Read the plot ranges
    mintime     = min(delays);
    maxtime     = max(delays);
    minabs      = min(rawsignal{j}(:));
    maxabs      = max(rawsignal{j}(:));
    zminmax     = round(max([abs(minabs) abs(maxabs)]),3);
    minwl       = min(cmprobe{j});
    maxwl       = max(cmprobe{j});
    Ncontours   = 40; % 40 contours by default is OK
    plotranges{j}  = [mintime maxtime minwl maxwl minabs maxabs Ncontours];
end
%% WRITE to dataStruct
% Write the main variables
dataStruct.delays      = delays;
dataStruct.cmprobe     = cmprobe;
dataStruct.rawsignal   = rawsignal;
dataStruct.noise       = noise;
dataStruct.plotranges  = plotranges;

% Calculate noise statistics
dataStruct.AvgNoise     = mean(noise(:));
dataStruct.MaxNoise     = max(noise(:));
dataStruct.SNR          = abs(round(zminmax/dataStruct.AvgNoise,3));

% Number of scans
   dataStruct.Nscans        = NaN;

for j=1:2
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

dataStruct.timescale = 'ps';