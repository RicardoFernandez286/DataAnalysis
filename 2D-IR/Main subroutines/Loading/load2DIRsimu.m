function  dataStruct = load2DIRsimu(dataStruct,varargin)

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
% Ricardo Fernandez-Teran / 18.08.2019 / v1.0a

%% DEBUG
% datafilename = 'R3_mTiO2_Re1213_MeOH_4PM_2D_154022';
% rootdir = 'D:\Data\20180418';
% dataStruct.rootdir=rootdir;

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
filename        = [rootdir filesep datafilename filesep 'spec_2D.dat'];

% Create progress bar and clear error string
if ShowWaitBar
    % dataStruct.WaitBar             = waitbar(0,'Loading data...');
    dataStruct.WBfigure                 = uifigure;
    dataStruct.WBfigure.Position(3:4)   = [405 175];
    dataStruct.WaitBar                  = uiprogressdlg(dataStruct.WBfigure,'Title','2D-IR data processing...','Message','Loading data...','Icon','info','ShowPercentage','on','Cancelable','on');
    dataStruct.ErrorText.String         = "";
end

if exist(filename,'file') ~= 0
Spec_2DIR_raw   = dlmread(filename);
delta_t2        = Spec_2DIR_raw(1,1);
Spec_2DIR_raw   = Spec_2DIR_raw(2:end,:);
Ndelays         = size(Spec_2DIR_raw,2) - 2;
t2delays        = delta_t2*(0:1:Ndelays-1)';
Nt1             = sqrt(size(Spec_2DIR_raw,1));
PROC_2D_DATA    = cell(Ndelays,1);
PumpAxis        = cell(Ndelays,1);
dummy_cell     = cell(Ndelays,1);

for i=1:Ndelays
    [W1,W3,Z]     = xyz2mat(Spec_2DIR_raw(:,[1:2 i+2]));
    PROC_2D_DATA{i,1}   = Z*1000;
    PumpAxis{i,1}       = W1;
    dummy_cell{i,1}     = 0;
    dummy_Onescell{i,1} = ones(Nt1,1); 
end

%% Write special variable
    dataStruct.isSimulation        = 1;
%% WRITE to dataStruct (Load)
    dataStruct.cmprobe       = W3;
    dataStruct.bins          = [];
    dataStruct.t2delays      = t2delays; % in ps!
    dataStruct.Ndelays       = Ndelays;
    dataStruct.Nspectra      = 1;
    dataStruct.Ndummies      = 1;
    dataStruct.Ndatastates   = 1;
    dataStruct.Nbins         = Nt1;
    dataStruct.Nslowmod      = 1;
    dataStruct.t1delays      = []; % in fs!
    dataStruct.Nscans        = 1;

%% WRITE to dataStruct (Process)
    dataStruct.ProbeAxis           = W3;
    dataStruct.freq_fit            = [];
    dataStruct.scattering_maxima   = [];
    dataStruct.PumpAxis            = PumpAxis;
    
    dataStruct.t1delays            = [];
    dataStruct.interferogram       = dummy_Onescell;
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
    dataStruct.apo_signal          = [];
    dataStruct.signal              = [];
    dataStruct.PROC_2D_DATA        = PROC_2D_DATA;
    dataStruct.SpecDiff            = 0;
    
% Clear 2DGC fit results
dataStruct.FitResults    = [];
dataStruct.t2_startFit   = [];
dataStruct.FitInput      = [];


else
    error('Not a valid 2D-IR dataset: empty folder or corrupt data')
end