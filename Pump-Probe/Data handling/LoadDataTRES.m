function dataStruct = LoadDataTRES(dataStruct)
%% READ from dataStruct
rootdir         = dataStruct.rootdir;
datafilename    = dataStruct.datafilename;

%% Check that all files exist and read them
fullName        = [rootdir filesep datafilename];

%% Prepare filenames
% Read the files if the directory is correctly populated
data            = csvread(fullName);

delays          = data(2:end,1);
cmprobe         = (data(1,2:end))';
rawsignal       = data(2:end,2:end);
noise           = zeros(size(rawsignal));
  
% Read the plot ranges
mintime     = min(delays);
maxtime     = max(delays);
minabs      = min(rawsignal(:));
maxabs      = max(rawsignal(:));
zminmax     = round(max([abs(minabs) abs(maxabs)]),3);
minwl       = min(cmprobe);
maxwl       = max(cmprobe);
Ncontours   = 40; % 40 contours by default is OK
plotranges  = [mintime maxtime minwl maxwl zminmax Ncontours minabs maxabs];

%% WRITE to dataStruct
% Write the main variables
dataStruct.delays      = delays;
dataStruct.cmprobe     = cmprobe;
dataStruct.rawsignal   = rawsignal;
dataStruct.noise       = noise;
dataStruct.plotranges  = plotranges;

% Enable all controls that are greyed out by default
dataStruct = EnableControls(dataStruct,'TRES');

% Set default values for certain controls
set(dataStruct.maxDeltaAbs_text,'String',[minabs maxabs]);
dataStruct.AddReplace1Dplot.Value = 1;
dataStruct.linlogtick.Value = 0;
dataStruct.SVD = 0;
percentwhites = 0;
dataStruct.percentwhites = percentwhites;
set(dataStruct.percent_whites,'String',num2str(percentwhites));

% Show noise statistics
AvgNoise = NaN;
MaxNoise = NaN;
SNR = NaN;
set(dataStruct.AvgNoise_text,'String',AvgNoise);
set(dataStruct.MaxNoise_text,'String',MaxNoise);
set(dataStruct.SNRnumber,'String',SNR);

% Set defaults for background subtraction subunit
dataStruct.BkgSubTick.Value = 0;
    % Show by default the background subtraction limits
    % from the first data point till just before time zero
    k = 1;
    t = dataStruct.delays(k);
    while t < -5
        k = k+1;
        t = dataStruct.delays(k+1);
    end
    dataStruct.j = 1;
    dataStruct.k = k;
    set(dataStruct.mintimeBkg,'String',num2str(dataStruct.delays(1)));
    set(dataStruct.maxtimeBkg,'String',num2str(dataStruct.delays(k)));
    % Do the background subtraction and change status in dataStruct.rawcorr
    dataStruct.bkg = mean(dataStruct.rawsignal(1:k,:));
    dataStruct.corrdata = dataStruct.rawsignal - dataStruct.bkg;
    dataStruct.rawcorr = 'RAW';
    
dataStruct.timescale = 'ns';