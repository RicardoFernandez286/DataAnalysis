function dataStruct = LoadData_nsTAUniGE(dataStruct)
%% READ from dataStruct
rootdir         = dataStruct.rootdir;
datafilename    = dataStruct.datafilename;

%% Check that all files exist and read them
fullName        = [rootdir filesep datafilename];

%% Decide whether to flip the sign of the signal
sign=-1;

%% Read Data
% Read the files if the directory is correctly populated
text            = fscanf(fopen(fullName),'%c');
head_n          = count(text,'%');

alldata         = readmatrix(fullName,'FileType','text','NumHeaderLines',head_n,'CommentStyle','%s');
delays          = unique(alldata(:,1))*1e9; % convert to ns

alldata         = alldata.*sign;
Ndelays         = length(delays);
Npixels         = size(alldata,1)./Ndelays;

rmvIdx          = delays >= 1e4;
delays(rmvIdx)  = [];
Ndelays         = length(delays);

%%%% remove NaN lines due to comments at end of file
% NaNidx          = sum(isnan(cmprobe));
if exist([rootdir filesep 'pix2lam.mat'],'file') ~= 0
    load([rootdir filesep 'pix2lam.mat'],'lam');
    cmprobe{1}  = lam;
else
    cmprobe{1}  = (1:Npixels)'; %#ok<*BDSCI>
end

%% Split the data by delays
tmpsignal   = zeros(Ndelays,Npixels);
tmpnoise    = zeros(Ndelays,Npixels);

for j=1:Ndelays
    idxPix          = (1:Npixels) + (j-1).*Npixels;
    cts(j)          = sign.*alldata(1+(j-1).*Npixels,6);
    tmpsignal(j,:)  = alldata(idxPix,3);
    tmpnoise(j,:)   = alldata(idxPix,5);
end

% figure;
% plot(delays,cts);

%% Remove points with zero counts
rmvIdx              = (cts == 0);
tmpsignal(rmvIdx,:) = [];
tmpnoise(rmvIdx,:)  = [];
delays(rmvIdx)      = [];
Ndelays             = length(delays);

% for p=1:Npixels
%     tmpsignal(:,p)  = smooth(delays,tmpsignal(:,p));
% end

% for j=1:Ndelays
%     tmpsignal(j,:)  = smooth(cmprobe{1},tmpsignal(j,:));
% end

%% Read single scan data
    Nscans      = NaN;
    noise{1}    = tmpnoise;
    rawsignal{1}= tmpsignal;

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
dataStruct.scandata    = [];
dataStruct.scanNoise   = [];

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

dataStruct.timescale = 'ns';
fclose('all');