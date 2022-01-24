function dataStruct = LoadDataIR5(dataStruct)
%% READ from dataStruct
% Get stuff from dataStruct
rootdir         = dataStruct.rootdir;
datafilename    = dataStruct.datafilename;


%% Read the file
data            = readmatrix([rootdir filesep datafilename]);
delays          = data(2:end,1)/1e3; % Delays originally in microseconds, convert to seconds
cmprobe         = data(1,2:end)';

rawsignal       = data(2:end,2:end);
Ndelays         = length(delays);

%% Interpolate and reduce number of datapoints
% DwF         = 400; % Downscaling factor
% Npoints     = round(Ndelays/DwF);
% NnegPts     = 100;
% InterpMethod= 'spline';
% 
% negDelays   = sort(-logspace(log10(min(-delays(delays<= 0))),log10(max(-delays)),NnegPts));
% posDelays   = logspace(log10(min(delays(delays>= 0))),log10(max(delays)),Npoints);
% newDelays   = [negDelays posDelays];
% rawsignal   = interp1(delays,rawsignal,newDelays,InterpMethod,'extrap');
% 
% delays      = newDelays';

%% Noise and ranges
noise           = zeros(size(rawsignal));

% Read the plot ranges
mintime         = min(delays);
maxtime         = max(delays);
minabs          = min(rawsignal(:));
maxabs          = max(rawsignal(:));
zminmax         = round(max([abs(minabs) abs(maxabs)]),3);
minwl           = min(cmprobe);
maxwl           = max(cmprobe);
Ncontours       = 40;
plotranges      = [mintime maxtime minwl maxwl minabs maxabs Ncontours];

%% WRITE to dataStruct
% Write the main variables
dataStruct.delays      = delays;
dataStruct.cmprobe     = {cmprobe};
dataStruct.rawsignal   = {rawsignal};
dataStruct.noise       = {noise};
dataStruct.plotranges  = {plotranges};

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
        dataStruct.bkg{1} = dataStruct.rawsignal{1}(1,:);
    else
        dataStruct.bkg{1} = mean(dataStruct.rawsignal{1}(Idx(1):Idx(2),:));
    end
dataStruct.corrdata{1} = dataStruct.rawsignal{1} - dataStruct.bkg{1};

dataStruct.timescale = 'ps';