function dataStruct = LoadData_nsTAUniGE(dataStruct)
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
delays          = unique(alldata(:,1))*1e9; % convert to ns
Ndelays         = length(delays);
Npixels         = size(alldata,1)./Ndelays;

removePix       = sort([1:25 Npixels-(0:1:10)]);

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
    cts(j)          = alldata(1+(j-1).*Npixels,6);
    tmpsignal(j,:)  = alldata(idxPix,3);
    tmpnoise(j,:)   = alldata(idxPix,5);
end
 
% signal = signal-signal(1,:);
% 
% maxPlt = 0.02.*max(signal(:));
% minPlt = -maxPlt;
% 
% ctrs = linspace(minPlt,maxPlt,40);
% contourf(cmprobe{1},delays,signal,ctrs,'EdgeColor','flat')
% cmap = darkb2r(minPlt,maxPlt,40,2);
% colormap(cmap);
% caxis([minPlt,maxPlt])
% ax=gca; ax.YScale='log';
% colorbar;
% 
% plot(delays,mean(signal(:,185:210),2),'o')
% ax.XScale='log';
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