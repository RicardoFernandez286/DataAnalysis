function  dataStruct = load2DIR_UoS(dataStruct,varargin)

% Description:  This function loads 2DIR data in the file format 
%               from the setup at the University of Sheffield
% Usage: dataStruct = load2DIR_UoS(dataStruct)
% Inputs:
%   dataStruct structure with fields:
%     datafilename
%     rootdir
%     
% Outputs:
%     cmprobe           (Double)
%     bins              (Double)
%     t2delays          (Double)
%     Ndelays           (Double)
%     Nspectra          (Double)
%     Ndatastates       (Double)
%     Nbins             (Double)
%     Nslowmod          (Double)
%     count             (Cell array)
%     probe             (Cell array)
%     reference         (Cell array)
%     interferogram     (Cell array)
%     signal			(Cell array)
%
% Ricardo Fernandez-Teran / 11.11.2021 / v2.5d

%% DEBUG
% rootdir = ''
% datafilename = ''
% ShowWaitBar = true;

%% HARDCODED Settings

%% READ from dataStruct
if isempty(varargin)
    ShowWaitBar = true;
else
    ShowWaitBar = false;
end

datafilename    = dataStruct.datafilename;
rootdir         = dataStruct.rootdir;
DetSz           = [96 96 32 32]; % Sizes of [probe1 probe2 ref1 ref2] in pixels

%% Load all the necessary files after checking that they exist
datadir     = [rootdir filesep datafilename];
datadir_fl  = dir(datadir);
datadir_fn  = {datadir_fl.name};

[~,fn,~]    = fileparts(datadir_fn{contains(datadir_fn(:),'.LG','IgnoreCase',1)});
filename    = [rootdir filesep datafilename filesep fn];

datadir_fn  = datadir_fn(~contains(datadir_fn(:),fn,'IgnoreCase',1))';

% Create progress bar and clear error string
if ShowWaitBar
    % dataStruct.WaitBar             = waitbar(0,'Loading data...');
    dataStruct.WBfigure                 = uifigure;
    dataStruct.WBfigure.Position(3:4)   = [405 175];
    dataStruct.WaitBar                  = uiprogressdlg(dataStruct.WBfigure,'Title','2D-IR data processing...','Message','Loading data...','Icon','info','ShowPercentage','on','Cancelable','on');
    dataStruct.ErrorText.String         = "";
end

% Read Spectrograph Information
text        = readlines([filename '.LG']);
g1_st       = strsplit(text{75},' '); % Det186 is 1st
w1_st       = strsplit(text{77},' ');
g2_st       = strsplit(text{67},' '); % Det185 is 2nd
w2_st       = strsplit(text{69},' ');

Gratings(1) = str2double(g1_st{end});
CWL(1)      = str2double(w1_st{end});
Gratings(2) = str2double(g2_st{end});
CWL(2)      = str2double(w2_st{end});

est_probe   = cell(1,2);

% Estimate probe axis, arbitrary calibration factors
for i=1:2
    switch Gratings(i)
        case 0 % 120 l/mm
            ppnm    = 0.168;
        case 1 % 100 l/mm
            ppnm    = 0.14;
        case 2 % 50 l/mm
            ppnm    = 0.068;
    end
    pix = flip(1:DetSz(i));
    est_probe{i} = 1e7./(pix./ppnm + CWL(i)+60 - DetSz(i)./2./ppnm);
end

% Read t2 delays
t2delays  = readmatrix([filename '.DT'],'FileType','delimitedtext');
Ndelays   = length(t2delays);

%% Build file list
fn_L        = datadir_fn(contains(datadir_fn(1:end),'.2D','IgnoreCase',1));
Ndone       = length(fn_L);
delName     = cell(Ndone,1);
delFN       = zeros(Ndone,1);

for i=1:Ndone
    delName{i}  = strsplit(fn_L{i},{'_D','_'});
    delFN(i)    = str2double(delName{i}{4});
end
[~,delIdx]  = sort(delFN);

Mcomplete   = (exist([filename '.2D'],'file') ~= 0);

if Mcomplete == 1 %if measurement is completed
    [t2delays,t2idx]  = sort(t2delays);
else
    [t2delays,~]    = sort(t2delays);
    t2delays        = t2delays(delIdx);
    Ndelays         = Ndone;
    [t2delays,t2idx]= sort(t2delays);
end


% If there is an w0 file in the current ROOTDIR, use it to set the rotating frame frequency.
% Otherwise, plot relative to w0
if exist([datadir filesep 'w0.csv'],'file') == 2
	rotframe = readmatrix([datadir filesep 'w0.csv']);
    w0  = rotframe(1);
    dt1 = rotframe(2);
else
    prompt = {'Enter rotating frame frequency (\omega_{0}, in cm^{-1}):','Enter t_{1} time step (in fs):'};
    dlgtitle = 'Rotating Frame Settings';
    definput = {'1700','22'};
    dims = [1 40];
    opts.Interpreter = 'tex';
    answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
    
    if isempty(answer)
        w0  = 0;
        dt1 = 1;
    else
        w0  = str2double(answer{1});
        dt1 = str2double(answer{2});
        try
            writematrix([w0;dt1],[datadir filesep 'w0.csv']);
        catch
%             warndlg('Error writing rotating frame info file','Error writing file');
            disp('Error writing rotating frame info file');
        end
    end
end

%% Read everything that has been done
signal_temp     = cell(Ndone,2);
t1delays        = cell(Ndone,2);
dummy_cell      = cell(Ndone,2);
dummy_Onescell  = cell(Ndone,2);
cmprobe         = cell(1,2);

if Ndone >= 1
for k=1:2 % Ndetectors
    idx = (1:DetSz(k)) + sum(DetSz(1:k-1));
    for i=1:Ndone
        rawdata             = readmatrix([datadir filesep fn_L{delIdx(i)}],'FileType','delimitedtext');
        if i==1
            Npixels             = size(rawdata,2);
            Nbins               = size(rawdata,1);
            bins                = (1:Nbins)'-1;
        end
        dummy_cell{i,k}     = 0;
        dummy_Onescell{i,k} = ones(Nbins,1);
%       signal_temp{i,k}    = -log10(rawdata((Nbins*(i-1)+1):(Nbins*i),idx)+1)*1000; % convert dT/T to mOD ->   -log10(dT/T + 1)*1000
        signal_temp{i,k}    = rawdata(:,idx); % data is in mOD - best for numerical accuracy when saving
        t1delays{i,k}       = [bins bins.*dt1];
    end
end

signal = signal_temp(t2idx,:);

%% WRITE to dataStruct (Load)
    dataStruct.Gratings      = Gratings;
    dataStruct.CWL           = CWL;
    dataStruct.est_probe     = est_probe; % Estimated probe from spectrograph info
    
    dataStruct.isSimulation  = 0;
    dataStruct.isShaper      = 1;
    dataStruct.cmprobe       = cmprobe;
    dataStruct.bins          = bins;
    dataStruct.t2delays      = t2delays; % Delays already in ps !!!
    dataStruct.Ndelays       = Ndelays;
    dataStruct.signal        = signal;
    dataStruct.Nspectra      = 2;
    dataStruct.Ndummies      = 1;
    dataStruct.Ndatastates   = 1;
    dataStruct.Nbins         = Nbins;
    dataStruct.Nslowmod      = 1;
    dataStruct.t1delays      = t1delays; % in fs!
    dataStruct.Nscans        = 1;
    dataStruct.datatype      = 'TimeFreq';
    dataStruct.interferogram = dummy_Onescell;
    dataStruct.w0            = w0;
    dataStruct.dt1           = dt1;
% Clear 2DGC fit results
dataStruct.FitResults    = [];
dataStruct.t2_startFit   = [];
dataStruct.FitInput      = [];

else
    error('Not a valid 2D-IR dataset: empty folder or corrupt data')
end
