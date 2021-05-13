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
% Ricardo Fernandez-Teran / 10.05.2021 / v0.9a

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

%% Load all the necessary files after checking that they exist
datadir     = [rootdir filesep datafilename];
datadir_fl  = dir(datadir);
datadir_fn  = {datadir_fl.name};

[~,fn,~]   = fileparts(datadir_fn{contains(datadir_fn(:),'.LG','IgnoreCase',1)});

filename    = [rootdir filesep datafilename filesep fn];

% Create progress bar and clear error string
if ShowWaitBar
    % dataStruct.WaitBar             = waitbar(0,'Loading data...');
    dataStruct.WBfigure                 = uifigure;
    dataStruct.WBfigure.Position(3:4)   = [405 175];
    dataStruct.WaitBar                  = uiprogressdlg(dataStruct.WBfigure,'Title','2D-IR data processing...','Message','Loading data...','Icon','info','ShowPercentage','on','Cancelable','on');
    dataStruct.ErrorText.String         = "";
end

if exist([filename '.2D'],'file') ~= 0
rawdata   = readmatrix([filename '.2D'],'FileType','delimitedtext');
t2delays  = readmatrix([filename '.DT'],'FileType','delimitedtext');

Ndelays   = length(t2delays);
Npixels   = size(rawdata,2);
Nbins     = round(size(rawdata,1)./Ndelays);
bins      = (1:Nbins)';

% Variable variables (!!)
dt1       = 14; % fs; may change
w0        = 1670; % need to figure this in a better way

% if exist([filename '.2D'],'file') ~= 0
%     cmprobe = 
% else
%     cmprobe   = (1:Npixels)';  % (?)
% end

signal          = cell(Ndelays,1);
t1delays        = cell(Ndelays,1);
dummy_cell      = cell(Ndelays,1);
dummy_Onescell  = cell(Ndelays,1);
for i=1:Ndelays
    dummy_cell{i,1}     = 0;
    dummy_Onescell{i,1} = ones(Nbins,1);
    signal{i,1}         = rawdata((Nbins*(i-1)+1):(Nbins*i),:);
    t1delays{i,1}       = bins.*dt1;
end

%% WRITE to dataStruct (Load)
    dataStruct.isSimulation  = 0;
    dataStruct.isShaper      = 1;
    dataStruct.cmprobe       = cmprobe;
    dataStruct.bins          = bins;
    dataStruct.t2delays      = t2delays; % Delays already in ps !!!
    dataStruct.Ndelays       = Ndelays;
    dataStruct.signal        = signal;
    dataStruct.Nspectra      = 1;
    dataStruct.Ndummies      = 1;
    dataStruct.Ndatastates   = 1;
    dataStruct.Nbins         = Nbins;
    dataStruct.Nslowmod      = 1;
    dataStruct.t1delays      = t1delays; % in fs!
    dataStruct.Nscans        = 1;
    dataStruct.datatype      = 'TimeFreq';
    dataStruct.interferogram = dummy_Onescell;

% Clear 2DGC fit results
dataStruct.FitResults    = [];
dataStruct.t2_startFit   = [];
dataStruct.FitInput      = [];

else
    error('Not a valid 2D-IR dataset: empty folder or corrupt data')
end