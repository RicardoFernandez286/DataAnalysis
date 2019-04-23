function dataStruct = LoadDataIR2(dataStruct,varargin)
%% READ from dataStruct
% Get stuff from dataStruct
rootdir         = dataStruct.rootdir;
datafilename    = dataStruct.datafilename;

if isempty(varargin)
    recalcscans = 0;
else
    recalcscans = varargin{1};
end
file_header     = [rootdir filesep datafilename filesep datafilename];


%% Prepare filenames
% Prepare the filenames for each file to read
delayfile       = [file_header '_delays.csv'];
wavenumberfile  = [file_header '_wavenumbers.csv'];
signalfile      = [file_header '_signal.csv'];
noisefile       = [file_header '_signal_noise.csv'];

% Read the files if the directory is correctly populated
existcondition = exist(delayfile, 'file') && exist(wavenumberfile, 'file') && exist(signalfile, 'file') && exist(noisefile,'file');
if existcondition == 1
% Read the files
delays          = csvread(delayfile);
cmprobe         = csvread(wavenumberfile);
if recalcscans == 1
    rawsignal       = dataStruct.rawsignal;
    noise           = dataStruct.noise;
else
    rawsignal       = csvread(signalfile);
    noise           = csvread(noisefile);
end
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
tempdir = [rootdir filesep datafilename filesep 'temp'];
if exist(tempdir,'dir') == 7
   filelist                 = dir([tempdir filesep '*signal*.csv']);
   dataStruct.Nscans        = floor(length(filelist)./2);
else
   dataStruct.Nscans        = NaN;
end

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
    
else
  assert(existcondition,'Selected directory is empty or files are corrupt!')
end
dataStruct.timescale = 'ns';