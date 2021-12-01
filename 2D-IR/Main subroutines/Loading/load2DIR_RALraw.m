% function  dataStruct = load2DIR_RALraw(dataStruct,varargin)

% Description: This function loads pre-processed 2DIR data from the ULTRA LIFEtime experiment at the
% Central Laser Facility of the Rutherford Appelton Laboratory.
%
% Usage: dataStruct = load2DIR_RALraw(dataStruct)
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
% Ricardo Fernandez-Teran / 13.05.2021 / v0.9a

%% DEBUG
rootdir         = 'D:\GoogleDrive\RESULTS\From James\2DIR DATA';
datafilename    = 'RAL LIFEtime - LT5 PTZCCC13C13NAP  2DIR 22 fs';
ShowWaitBar     = true;

%% HARDCODED Settings

%% READ from dataStruct
% if isempty(varargin)
    ShowWaitBar = true;
% else
%     ShowWaitBar = false;
% end
% 
% datafilename    = dataStruct.datafilename;
% rootdir         = dataStruct.rootdir;

if 1==1
%% Load all the necessary files after checking that they exist
datadir     = [rootdir filesep datafilename];
datadir_fl  = dir(datadir);
datadir_fn  = {datadir_fl.name};

[~,fn,~]    = fileparts({datadir_fn{contains(datadir_fn(:),'.csv','IgnoreCase',1)}});
fn          = fn(contains(fn,'Avg Diff','ignorecase',1));

alldata     = dlmread([datadir filesep fn{1} '.csv'],',',0,1);

%%
rawdata     = alldata(2:end,:);

delayarray  = alldata(1,:);
t1ranges    = [0 find(diff(delayarray) < 0) size(rawdata,2)];
t2idx       = [1 find(diff(delayarray) < 0)];

dt1         = 22; % fs
Nphases     = 4;

Nbins       = unique(diff(t1ranges))/Nphases;
if length(Nbins) > 1
    error('Irregular number of t1 delays... Not yet implemented! Aborting.');
elseif mod(Nbins,1) ~= 0
    error('Inconsistent number of bins for phase cycling, please check!');
end

t2delays    = delayarray(t2idx);
bins        = (1:Nbins)';
Ndelays     = length(t2delays);
Npixels     = size(alldata,1)-1;

% clear alldata

signal4D    = zeros(Npixels,Ndelays,Nbins,Nphases);

for i=1:Ndelays
    binIdx_del      = (t1ranges(i)+1):t1ranges(i+1);
    binIdx_phC      = [0 (1:4)*Nbins];
    for j=1:4
        signal4D(:,i,:,j) = rawdata(:,binIdx_del((binIdx_phC(j)+1):binIdx_phC(j+1))); 
    end
end

%% Calculate Phase Cycling
% signal = + 1 - 2 + 3 - 4

k=1;
signal = cell(Ndelays,1);

phase_sig   = [+1 -1 +1 -1];
sig     = zeros(Npixels,Nbins,Nphases);
% scatt   = zeros(Npixels,Nbins,Nphases);

for m=1:Ndelays
    for p=1:Nphases
        sig(:,:,p)     = squeeze(signal4D(:,m,:,p)).*phase_sig(p);
%         scatt(:,:,p)   = squeeze(signal4D(:,m,:,p)).*phase_scatt(p);
    end
    signal{m,k} = squeeze(sum(sig,3));
%     scattering{m,k}  = squeeze(sum(scatt,3));
end

%%
% Create progress bar and clear error string
if ShowWaitBar == 1
    progress = 0;
    % dataStruct.WaitBar             = waitbar(0,'Loading data...');
    dataStruct.WBfigure                 = uifigure;
    dataStruct.WBfigure.Position(3:4)   = [405 175];
    dataStruct.WaitBar                  = uiprogressdlg(dataStruct.WBfigure,'Title','2D-IR data processing...','Message','Loading data...','Icon','info','ShowPercentage','on','Cancelable','on');
    dataStruct.ErrorText.String         = "";
end


Npixels   = length(cmprobe);

lowpart   = 1:(Npixels/2);
highpart  = (Npixels/2:Npixels-1)+1;

ProbeCut  = 'other';

switch ProbeCut
    case 'low'
        startlow    = 1;
        endlow      = length(cmprobe(cmprobe(lowpart) < min(cmprobe(highpart))));
        starthigh   = Npixels/2+1;
        endhigh     = Npixels;
    case 'high'
        startlow    = 1;
        endlow      = Npixels/2;
        starthigh   = Npixels/2 + find(cmprobe(highpart) > max(cmprobe(lowpart)),1);
        endhigh     = Npixels;
    case 'other'
        startlow    = 1;
        endlow      = Npixels/2 - 20;
        starthigh   = Npixels/2 + find(cmprobe(highpart) > max(cmprobe(startlow:endlow)),1);
        endhigh     = Npixels;
    case 'none'
        startlow    = 1;
        endlow      = Npixels/2;
        starthigh   = Npixels/2 +1;
        endhigh     = Npixels;
end
idProbe = [startlow:endlow starthigh:endhigh];
Npixels = length(idProbe);
cmprobe = cmprobe(idProbe);
% idProbe   = 1:Npixels;

signal          = cell(Ndelays,1);
t1delays        = cell(Ndelays,1);
dummy_cell      = cell(Ndelays,1);
dummy_Onescell  = cell(Ndelays,1);

for i=1:Ndelays
    if ShowWaitBar == 1
        % Update the Wait Bar
        progress = progress + 1;
        dataStruct.WaitBar.Value    = (progress/(Ndelays));
        dataStruct.WaitBar.Message  = ['Processing data... (' num2str(progress) ' of ' num2str(Ndelays) ')'];
        if dataStruct.WaitBar.CancelRequested
            delete(dataStruct.WBfigure);
            error('User aborted loading the data!');
        end
        drawnow;
    end
    dummy_cell{i,1}     = 0;
    dummy_Onescell{i,1} = ones(Nbins,1);
    signal{i,1}         = zeros(Nbins,Npixels);
    t1delays{i,1}       = bins*dt1;
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

    return
%% WRITE to dataStruct (Process 2D-IR)
    dataStruct.ProbeAxis           = cmprobe;
    dataStruct.freq_fit            = [];
    dataStruct.scattering_maxima   = [];
    dataStruct.PumpAxis            = PumpAxis;
    
    dataStruct.t1delays            = t1delays;
    
    dataStruct.binzero             = dummy_cell;
    dataStruct.binspecmax          = ones(Ndelays,1);
    dataStruct.apodize_function    = dummy_Onescell;
    dataStruct.FFT_ZPsig           = [];
    dataStruct.phased_FFTZPsig     = dummy_Onescell;
	dataStruct.phased_FFTZPint     = dummy_Onescell;
    dataStruct.fittedPhase         = dummy_Onescell;
    dataStruct.phasepoints         = dummy_Onescell;
    dataStruct.ZP_phase            = dummy_cell;
    dataStruct.phase_coeff         = dummy_cell;
    dataStruct.apo_interferogram   = [];
    dataStruct.apo_signal          = signal;
    dataStruct.PROC_2D_DATA        = PROC_2D_DATA;
    dataStruct.SpecDiff            = 0;
    
% Clear 2DGC fit results
dataStruct.FitResults    = [];
dataStruct.t2_startFit   = [];
dataStruct.FitInput      = [];

else
    error('Not a valid 2D-IR dataset: empty folder or corrupt data')
end