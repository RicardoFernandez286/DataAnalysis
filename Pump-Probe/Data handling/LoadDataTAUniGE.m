function dataStruct = LoadDataTAUniGE(dataStruct)
%% READ from dataStruct
rootdir         = dataStruct.rootdir;
datafilename    = dataStruct.datafilename;

%% Check that all files exist and read them
fullName        = [rootdir filesep datafilename];

%% Read Data
% Read the files if the directory is correctly populated
text            = fscanf(fopen(fullName),'%c');
head_n          = count(text,'%');

alldata         = readmatrix(fullName,'FileType','text','NumHeaderLines',head_n,'CommentStyle','%s');

Npixels         = round(size(alldata,2)/2)-1;

removePix       = sort([1:25 Npixels-(0:1:10)]);

% cmprobe         = 
delays          = unique((alldata(2:end,1)))*1e12; % convert to ps
Ndelays         = length(delays);

%%%% remove NaN lines due to comments at end of file
% NaNidx          = sum(isnan(cmprobe));
if exist([rootdir filesep 'pix2lam.mat'],'file') ~= 0
    load([rootdir filesep 'pix2lam.mat'],'lam');
    cmprobe     = lam;
else
    cmprobe     = (1:Npixels)';
end

%% Read single scan data
if dataStruct.chirpCorr == 0
    Nscans          = round(size(alldata,1)/length(delays));
    rawsignal       = zeros(Ndelays,Npixels);
else
    % Data is chirp-corrected, do NOT reload from files
    rawsignal       = dataStruct.rawsignal;
    cmprobe         = dataStruct.cmprobe;
    delays          = dataStruct.delays;
    Nscans          = 0;
end

if Nscans >= 1
    scandata            = zeros([size(rawsignal) Nscans]);
    scan_err            = zeros([size(rawsignal) Nscans]);
    for s=1:Nscans
        scandata(:,:,s) = alldata((s-1)*Ndelays+(1:Ndelays),3:2:(2*Npixels+2))./1e3; % s-th set of delays (rows), odd columns
        scan_err(:,:,s) = alldata((s-1)*Ndelays+(1:Ndelays),4:2:(2*Npixels+2))./1e3; % s-th set of delays (rows), even columns
    end
    noise                   = mean(scan_err,3); % in mOD
    rawsignal               = mean(scandata,3); % in mOD
    scandata(:,removePix,:) = [];
    scan_err(:,removePix,:) = [];
    rawsignal(:,removePix)  = [];
    noise(:,removePix)      = [];
    cmprobe(removePix)      = [];
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
dataStruct.scandata    = scandata;
dataStruct.scanNoise   = scan_err;

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
fclose('all');