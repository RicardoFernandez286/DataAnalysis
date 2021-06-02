function dataStruct = LoadDataIR_UoS(dataStruct,varargin)
%% READ from dataStruct
% Get stuff from dataStruct
rootdir         = dataStruct.rootdir;
datafilename    = dataStruct.datafilename;

% Get files
datadir     = [rootdir filesep datafilename];
datadir_fl  = dir(datadir);
datadir_fn  = {datadir_fl.name};

[~,fn,~]   = fileparts(datadir_fn{contains(datadir_fn(:),'.LG','IgnoreCase',1)});

filename    = [rootdir filesep datafilename filesep fn];

%% Read the file
data            = readmatrix([filename '.2D'],'FileType','delimited');
delays          = readmatrix([filename '.DT'],'FileType','delimited'); % Delays originally in microseconds, convert to seconds
Npixels         = 96;
cmprobe         = 1:Npixels;

rawsignal       = data;
Ndelays         = length(delays);

%% Noise and ranges
noise           = zeros(size(rawsignal));

% Read the plot ranges
mintime     = min(delays);
maxtime     = max(delays);
minabs      = min(rawsignal(:));
maxabs      = max(rawsignal(:));
zminmax     = round(max([abs(minabs) abs(maxabs)]),3);
minwl       = min(cmprobe);
maxwl       = max(cmprobe);
plotranges  = [mintime maxtime minwl maxwl minabs maxabs];

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

% Background Subtraction
if dataStruct.recalcBkg == 0
    % Subtract the first negative delay from all the dataset by default - otherwise take the inputs
    dataStruct.mintimeBkg = dataStruct.delays(1);
    dataStruct.maxtimeBkg = dataStruct.delays(1);
end
    Idx = findClosestId2Val(dataStruct.delays,[dataStruct.mintimeBkg dataStruct.maxtimeBkg]);
    % Do the background subtraction and change status in handles.rawcorr
    if Idx(2)==1
        dataStruct.bkg = dataStruct.rawsignal(1,:);
    else
        dataStruct.bkg = mean(dataStruct.rawsignal(Idx(1):Idx(2),:));
    end
dataStruct.corrdata = dataStruct.rawsignal - dataStruct.bkg;

dataStruct.timescale = 's';