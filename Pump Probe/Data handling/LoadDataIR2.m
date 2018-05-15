function handles = LoadDataIR2(handles)
%% READ from handles
% Get rootdir and datafilename from handles
rootdir         = char(handles.CurrDir.String);
datafilename    = char(handles.datafilename);

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
rawsignal       = csvread(signalfile);
noise           = csvread(noisefile);
  
% Read the plot ranges
mintime = min(delays);
maxtime = max(delays);
minabs = min(rawsignal(:));
maxabs = max(rawsignal(:));
zminmax = round(max([abs(minabs) abs(maxabs)]),3);
minwl = min(cmprobe);
maxwl = max(cmprobe);
Ncontours = 40; % 40 contours by default is OK
plotranges = [mintime maxtime minwl maxwl zminmax Ncontours minabs maxabs];

%% WRITE to handles
% Write the main variables
handles.delays = delays;
handles.cmprobe = cmprobe;
handles.rawsignal = rawsignal;
handles.noise = noise;
handles.plotranges = plotranges;

% Enable all controls that are greyed out by default
handles = EnableControls(handles,'IRLab2');

% Set default values for certain controls
set(handles.maxDeltaAbs_text,'String',[minabs maxabs]);
handles.AddReplace1Dplot.Value = 1;
handles.SVD = 0;
percentwhites = 5;
handles.percentwhites = percentwhites;
set(handles.percent_whites,'String',num2str(percentwhites));

% Show noise statistics
AvgNoise = round(mean(mean(noise)),3);
MaxNoise = round(max(max(noise)),3);
SNR = round(zminmax/AvgNoise,3);
set(handles.AvgNoise_text,'String',AvgNoise);
set(handles.MaxNoise_text,'String',MaxNoise);
set(handles.SNRnumber,'String',SNR);

% Set defaults for background subtraction subunit
handles.BkgSubTick.Value = 1;
    % Show by default the background subtraction limits
    % from the first data point till just before time zero
    k = 1;
    t = handles.delays(k);
    while t < 0
        k = k+1;
        t = handles.delays(k+1);
    end
    handles.j = 1;
    handles.k = k;
    set(handles.mintimeBkg,'String',num2str(handles.delays(1)));
    set(handles.maxtimeBkg,'String',num2str(handles.delays(k)));
    % Do the background subtraction and change status in handles.rawcorr
    handles.bkg = mean(handles.rawsignal(1:k,:));
    handles.corrdata = handles.rawsignal - handles.bkg;
    handles.rawcorr = 'CORRECTED';
else
  assert(existcondition,'Selected directory is empty or files are corrupt!')
end
handles.timescale = 'ns';