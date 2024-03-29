function dataStruct = LoadDataTAUniGE(dataStruct,recalcScans)
%% READ from dataStruct
rootdir         = dataStruct.rootdir;
datafilename    = dataStruct.datafilename;

%% Check that all files exist and read them
fullName        = [rootdir filesep datafilename];

%% Decide whether to flip the sign of the signal
sign=+1;

%% Read Data
% Read the files if the directory is correctly populated
text            = fscanf(fopen(fullName),'%c');
head_n          = count(text,'%') - count(text,'%%');
alldata         = readmatrix(fullName,'FileType','text','NumHeaderLines',head_n,'CommentStyle','%s');

Npixels         = round(size(alldata,2)/2)-1;
removePix       = sort([1:25 Npixels-(0:1:10)]);

alldelays       = alldata(1:end,1)*1e12; % convert to ps
[delays,~,dID]  = unique(alldelays); 
Ndelays         = length(delays);

alldata         = alldata.*sign;

%%%% remove NaN lines due to comments at end of file [currently not needed]
% NaNidx          = sum(isnan(cmprobe));

% Read the pixel to lambda calibration if it exists
if exist([rootdir filesep 'pix2lam.mat'],'file') ~= 0
    load([rootdir filesep 'pix2lam.mat'],'lam');
    cmprobe{1}  = lam;
else
    cmprobe{1}  = (1:Npixels)';
end

if cmprobe{1}(1) >= cmprobe{1}(end)
    removePix =  sort([1:50 Npixels-(0:1:95)]);
end

%% Read single scan data
if dataStruct.chirpCorr == 0
    Nscans          = floor(size(alldata,1)/Ndelays);
    rawsignal{1}    = zeros(Ndelays,Npixels);
    
    % remove NaN lines due to comments at end of file
    NaNidx          = sum(isnan(cmprobe{1}));
    cmprobe{1}      = cmprobe{1}(1:end-NaNidx);
    rawsignal{1}    = rawsignal{1}(:,1:end-NaNidx).*1e3; % convert to mOD
else
    % Data is chirp-corrected or scans recalculation, do NOT reload from files
    rawsignal       = dataStruct.rawsignal;
    noise           = dataStruct.noise;
    cmprobe         = dataStruct.cmprobe;
    delays          = dataStruct.delays;
    Nscans          = 0;
end

if recalcScans == 1
    rawsignal       = dataStruct.rawsignal;
    noise           = dataStruct.noise;
    cmprobe         = dataStruct.cmprobe;
    delays          = dataStruct.delays;
    Nscans          = 0;
else
    % remove NaN lines due to comments at end of file
    NaNidx          = sum(isnan(cmprobe{1}));
    cmprobe{1}      = cmprobe{1}(1:end-NaNidx);
    rawsignal{1}    = rawsignal{1}(:,1:end-NaNidx)*1e3; % convert to mOD
    rawsignal{1}    = fillmissing(rawsignal{1}, 'linear');
end
  
if Nscans >= 1
    scandata{1}         = zeros([size(rawsignal{1}) Nscans]);
    scan_err{1}         = zeros([size(rawsignal{1}) Nscans]);
    for s=1:Nscans
        IDs                 = (s-1)*Ndelays + (1:Ndelays);
        tmpScanData         = alldata(IDs,3:2:(2*Npixels+2))./1e3;
        smpScanError        = alldata(IDs,4:2:(2*Npixels+2))./1e3;
        scandata{1}(:,:,s)  = tmpScanData(dID(IDs),:);  % s-th set of delays (rows), odd columns     - indexed by dID [to sort into correct delay]
        scan_err{1}(:,:,s)  = smpScanError(dID(IDs),:); % s-th set of delays (rows), even columns   - indexed by dID [to sort into correct delay]
    end
    noise{1}                    = mean(scan_err{1},3); % in mOD
    rawsignal{1}                = mean(scandata{1},3); % in mOD
%     rawsignal{1}                = fillmissing(rawsignal{1}, 'linear');
    scandata{1}(:,removePix,:)  = [];
    scan_err{1}(:,removePix,:)  = [];
    rawsignal{1}(:,removePix)   = [];
    noise{1}(:,removePix)       = [];
    cmprobe{1}(removePix)       = [];
else
    Nscans      = NaN;
    noise{1}    = zeros(size(rawsignal{1}));
    scandata    = [];
    scan_err    = [];
end

if recalcScans ~= 1
    % Flip the data and wavelength axes if it's not increasing (e.g. NIR data)
    if cmprobe{1}(1) >= cmprobe{1}(end)
        cmprobe{1}   = flipud(cmprobe{1});
        rawsignal{1} = fliplr(rawsignal{1});
        noise{1}     = fliplr(rawsignal{1});
    end
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