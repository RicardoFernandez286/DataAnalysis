function  dataStruct = load2DIR_P2DAT(dataStruct,varargin)

% Description: This function loads pre-processed 2DIR data that has been saved in the P2DAT format (XYZZZZZZ)
%
% Usage: dataStruct = load2DIR_P2DAT(dataStruct)
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
% Ricardo Fernandez-Teran / 05.09.2022 / v1.0a

%% HARDCODED Settings
dt1 = 1; % Since the data is already in frequency domain

%% READ from dataStruct
if isempty(varargin)
    ShowWaitBar = true;
else
    ShowWaitBar = false;
end

datafilename    = dataStruct.datafilename;
rootdir         = dataStruct.rootdir;

%% Load Data
dataXYZZ = readmatrix([rootdir filesep datafilename],'FileType','text');
Ndelays  = size(dataXYZZ,2) - 2; % 1st two columns contain the pump and probe axes

[pump,probe,~] = xyz2mat(dataXYZZ(:,1:3));

w1      = pump(2:end);
w3      = probe(2:end);

Nbins   = length(w1);
Npix    = length(w3);

if Ndelays > 0
bins            = (1:Nbins)';

signal          = cell(Ndelays,2);
t1delays        = cell(Ndelays,2);
dummy_cell      = cell(Ndelays,1);
dummy_Onescell  = cell(Ndelays,1);
PumpAxis        = cell(Ndelays,1);
PROC_2D_DATA    = cell(Ndelays,1);

progress = 0;

% Create progress bar and clear error string
if ShowWaitBar == 1
    progress = 0;
    % dataStruct.WaitBar             = waitbar(0,'Loading data...');
    dataStruct.WBfigure                 = uifigure;
    dataStruct.WBfigure.Position(3:4)   = [405 175];
    dataStruct.WaitBar                  = uiprogressdlg(dataStruct.WBfigure,'Title','2D-IR data processing...','Message','Loading data...','Icon','info','ShowPercentage','on','Cancelable','on');
    dataStruct.ErrorText.String         = "";
end

for m=1:Ndelays 
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
    
    k=1; % only ONE detector since we're dealing with P2DAT files
    dummy_cell{m,k}     = 0;
    dummy_Onescell{m,k} = ones(Nbins,1);
    binzero{m,k}        = 1; % Because we are using the shaper
    binspecmax(m,k)     = 1; %#ok<*AGROW> 
        
    signal{m,k}         = zeros(Nbins,Npix);
    t1delays{m,k}       = [bins bins.*dt1];
    PumpAxis{m,k}       = w1;
    PROC_2D_DATA{m,k}   = reshape(dataXYZZ(2:end,m+2),[Nbins Npix]);
end
    
    ProbeAxis = {w3'};
    t2delays  = dataXYZZ(1,3:end)';
%% WRITE to dataStruct (Load)
    dataStruct.isSimulation  = -1; % Hack to deal with P2DAT files!
    dataStruct.isShaper      = 1;
    dataStruct.cmprobe       = ProbeAxis;
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
    dataStruct.datatype      = 'FreqFreq';
    dataStruct.interferogram = dummy_Onescell;

%% WRITE to dataStruct (Process 2D-IR)
    dataStruct.ProbeAxis           = ProbeAxis;
    dataStruct.freq_fit            = dummy_Onescell;
    dataStruct.scattering_maxima   = dummy_Onescell;
    dataStruct.PumpAxis            = PumpAxis;
    
    dataStruct.t1delays            = t1delays;
    
    dataStruct.binzero             = binzero;
    dataStruct.binspecmax          = binspecmax;
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