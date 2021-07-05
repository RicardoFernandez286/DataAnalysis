function dataStruct = LoadDataHelios(dataStruct)
%% READ from dataStruct
rootdir         = dataStruct.rootdir;
datafilename    = dataStruct.datafilename;

%% Check that all files exist and read them
fullName        = [rootdir filesep datafilename];

%% Read Data
% Read the files if the directory is correctly populated
data            = readmatrix(fullName,'CommentStyle','%s')';

delays          = data(2:end,1);
cmprobe         = (data(1,2:end))';
rawsignal       = data(2:end,2:end);

% remove NaN lines due to comments at end of file
NaNidx          = sum(isnan(cmprobe));
cmprobe         = cmprobe(1:end-NaNidx);
rawsignal       = rawsignal(:,1:end-NaNidx)*1e3; % convert to mOD
rawsignal       = fillmissing(rawsignal, 'linear');

%% Read single scan data
% dataname        = strsplit(datafilename,' ');
% [~,shortName]   = fileparts(dataname{1});
[~,shortName]   = fileparts(datafilename);
filelist        = dir(rootdir);
scanNames       = filelist(contains({filelist.name},[shortName '_scan'],'ignorecase',1) & contains({filelist.name},'csv','ignorecase',1));
scanNames       = {scanNames.name}';
Nscans          = length(scanNames);

if Nscans >= 1
    scandata            = zeros([size(rawsignal) Nscans]);
    for s=1:Nscans
        tempdata        = readmatrix([rootdir filesep scanNames{s}],'CommentStyle','%s')';
        NaNidx          = sum(isnan(tempdata(1,:)));
        tempdelays      = tempdata(2:end,1);
        tempWL          = tempdata(1,2:end);
        tempWL          = tempWL(1:end-NaNidx);
        tempscandata    = tempdata(2:end,2:(end-NaNidx));
        scandata(:,:,s) = fillmissing(tempscandata, 'linear');
    end
    noise               = std(scandata,0,3);
else
    Nscans = NaN;
    noise  = zeros(size(rawsignal));
end

%% Read the plot ranges
mintime     = min(delays);
maxtime     = max(delays);
minabs      = min(rawsignal(:));
maxabs      = max(rawsignal(:));
zminmax     = round(max([abs(minabs) abs(maxabs)]),3);
minwl       = min(cmprobe);
maxwl       = max(cmprobe);
Ncontours   = 40; % 40 contours by default is OK
plotranges  = [mintime maxtime minwl maxwl minabs maxabs Ncontours];

%% WRITE to dataStruct
% Write the main variables
dataStruct.delays      = delays;
dataStruct.cmprobe     = cmprobe;
dataStruct.rawsignal   = rawsignal;
dataStruct.noise       = noise;
dataStruct.plotranges  = plotranges;

% Calculate noise statistics
dataStruct.AvgNoise    = mean(noise(:),'omitnan');
dataStruct.MaxNoise    = max(noise(:));
dataStruct.SNR         = abs(round(zminmax/dataStruct.AvgNoise,3));

% Number of scans
dataStruct.Nscans      = Nscans;

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

dataStruct.timescale = 'ps';