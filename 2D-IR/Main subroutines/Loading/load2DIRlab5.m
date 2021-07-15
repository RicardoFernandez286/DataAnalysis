function  dataStruct = load2DIRlab5(dataStruct,varargin)

% Description: This function loads 2DIR data in the file format from the VET simulations
% Usage: dataStruct = load2DIRsimu(dataStruct)
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
% Ricardo Fernandez-Teran / 02.02.2021 / v1.0a

%% DEBUG
% datafilename = 'data2D_avg';
% rootdir = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\RESULTS\2D-IR\Lab5';
% ShowWaitBar = false;

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
filename        = [rootdir filesep datafilename '.dat'];

% Create progress bar and clear error string
if ShowWaitBar
    % dataStruct.WaitBar             = waitbar(0,'Loading data...');
    dataStruct.WBfigure                 = uifigure;
    dataStruct.WBfigure.Position(3:4)   = [405 175];
    dataStruct.WaitBar                  = uiprogressdlg(dataStruct.WBfigure,'Title','2D-IR data processing...','Message','Loading data...','Icon','info','ShowPercentage','on','Cancelable','on');
    dataStruct.ErrorText.String         = "";
end

if exist(filename,'file') ~= 0
rawdata   = readmatrix(filename);

cmprobe   = rawdata(1,2:end);

Ndelays   = unique(sum(rawdata(:,2:end)==cmprobe));
Nbins     = round((size(rawdata,1)-Ndelays)./Ndelays);
bins      = 1:Nbins;

newt2idx  = logical(sum(rawdata(:,2:end)==cmprobe,2)./length(cmprobe));

t2delays  = rawdata(newt2idx,1);

startidx        = find(newt2idx);
signal          = cell(Ndelays,1);
t1delays        = cell(Ndelays,1);
dummy_cell      = cell(Ndelays,1);
dummy_Onescell  = cell(Ndelays,1);
for i=1:Ndelays
    dummy_cell{i,1}     = 0;
    dummy_Onescell{i,1} = ones(Nbins,1);
    signal{i,1}         = rawdata(startidx(i)+bins,2:end)./40;
    t1delays{i,1}       = [(bins-1)' rawdata(startidx(i)+bins,1)];
end

%% WRITE to dataStruct (Load)
    dataStruct.isSimulation  = 0;
    dataStruct.isShaper      = 1;
    dataStruct.cmprobe       = cmprobe';
    dataStruct.bins          = bins;
    dataStruct.t2delays      = t2delays./1000; % in ps!
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
else
    error('Not a valid 2D-IR dataset: empty folder or corrupt data')
end