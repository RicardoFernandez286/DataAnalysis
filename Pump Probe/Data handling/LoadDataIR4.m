function handles = LoadDataIR4(handles)

% %% Debug
% ID = 1;

%% READ from handles
ID = handles.DatafoldersList.Value;
scanindices = handles.scanindices;
scandata    = handles.scandata;
alldelays   = handles.alldelays;
cmprobe     = handles.cmprobe;
SampleName  = handles.SampleName;
DelayZero   = handles.DelayZero;

%% Calculate the signal average and the standard deviation ("noise") for the selected data file

if size(alldelays,2) > 1
    delays = mean(alldelays(:,scanindices(ID,1):scanindices(ID,2)),2);
    delays = delays - DelayZero(ID);
end
rawsignal = mean(scandata(:,:,scanindices(ID,1):scanindices(ID,2)),3);
noise = std(scandata(:,:,scanindices(ID,1):scanindices(ID,2)),0,3);

%% Read the plot ranges
mintime = min(delays);
maxtime = max(delays);
minabs = min(rawsignal(:));
maxabs = max(rawsignal(:));
zminmax = round(max([abs(minabs) abs(maxabs)]),3);
minwl = min(cmprobe);
maxwl = max(cmprobe);
Ncontours = 20; % 20 contours by default is OK
plotranges = [mintime maxtime minwl maxwl zminmax Ncontours minabs maxabs];

% Write to handles
handles.delays = delays;
handles.rawsignal = rawsignal;
handles.noise = noise;
handles.plotranges = plotranges;
handles.datafilename = SampleName{ID};

% Enable all controls that are greyed out by default
handles = EnableControls(handles,'IRLab2');

% Set default values for certain controls
set(handles.maxDeltaAbs_text,'String',zminmax);
handles.AddReplace1Dplot.Value = 1;
handles.SVD = 0;
percentwhites = 2;
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
handles.timescale = 'ps';