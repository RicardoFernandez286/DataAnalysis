function dataStruct = LoadPDAT(dataStruct)
%% READ from dataStruct
rootdir         = dataStruct.rootdir;
datafilename    = dataStruct.datafilename;

%% Check that all files exist and read them
fullName        = [rootdir filesep datafilename];

%% Hardcoded settings
convertTimescale = 0; % Will convert ps data to ns

%% Read Data
% Read the files if the directory is correctly populated
fid             = fopen(fullName);
    units = fgetl(fid);
fclose(fid);
units           = strsplit(units,'*');

timescale       = units{1};
Xunits          = units{2};
switch Xunits
    case 'nm'
        probeunits = 'Wavelength';
    case 'cm^{-1}'
        probeunits = 'Wavenumbers';
    otherwise
        probeunits = 'Pixels';
end

data            = readmatrix(fullName,'FileType','text','CommentStyle','%s');

delays          = data(2:end,1);
cmprobe         = (data(1,2:end))';
rawsignal       = data(2:end,2:end);

% remove NaN lines due to comments at end of file
NaNidx          = sum(isnan(cmprobe));
cmprobe         = cmprobe(1:end-NaNidx);
rawsignal       = rawsignal(:,1:end-NaNidx); % convert to mOD
rawsignal       = fillmissing(rawsignal, 'linear');

Nscans = NaN;
noise  = zeros(size(rawsignal));

% Convert ps to ns if requested
if convertTimescale == 1
    switch timescale
        case 'ps'
            timescale = 'ns';
            delays = delays./1000;
            warning('Delays have been converted to ns!')
    end
end
%% Read the plot ranges
mintime     = min(delays); %#ok<*NASGU> 
maxtime     = max(delays);
minabs      = min(rawsignal(:));
maxabs      = max(rawsignal(:));
zminmax     = round(max([abs(minabs) abs(maxabs)]),3);
minwl       = min(cmprobe);
maxwl       = max(cmprobe);
Ncontours   = 40; % 40 contours by default is OK
plotranges  = [0.01 maxtime minwl maxwl minabs maxabs Ncontours];

%% WRITE to dataStruct
% Write the main variables
dataStruct.delays      = delays;
dataStruct.cmprobe     = {cmprobe};
dataStruct.rawsignal   = {rawsignal};
dataStruct.noise       = {noise};
dataStruct.plotranges  = {plotranges};

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
    dataStruct.bkg      = {};
    dataStruct.corrdata = {};
    Idx = findClosestId2Val(dataStruct.delays,[dataStruct.mintimeBkg dataStruct.maxtimeBkg]);
    % Do the background subtraction and change status in handles.rawcorr
    if Idx(2)==1
        dataStruct.bkg{1}  = dataStruct.rawsignal{1}(1,:);
    else
        dataStruct.bkg{1}  = mean(dataStruct.rawsignal{1}(Idx(1):Idx(2),:));
    end
dataStruct.corrdata{1}     = dataStruct.rawsignal{1} - dataStruct.bkg{1};

dataStruct.timescale    = timescale;
dataStruct.Xunits       = Xunits;
dataStruct.probeunits   = probeunits;