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

%% Read 1D spectra
    LinSpec_FD  = readmatrix([rootdir filesep datafilename filesep 'spec_lin.dat'],'Delimiter','  ');

%% Read surface trajectories
    % Read the first line
    trajectInfo     = rmmissing(readmatrix([rootdir filesep datafilename filesep 'trajectory.dat'],'Range',[1 1 1 6]));
    switch length(trajectInfo)
        case 6
            version = 3;
        case 4
            version = 2;
        case 3
            version = 1;
    end
    % Read the rest
    trajectData_2D  = readmatrix([rootdir filesep datafilename filesep 'trajectory.dat'],'NumHeaderLines',1);
    if version > 2
        Nmolecules  = groupcounts(trajectData_2D(:,8));
        Nsamples    = max(trajectData_2D(:,8));
    else
        Nmolecules  = trajectInfo(1);
        Nsamples    = size(trajectData_2D,1)/Nmolecules;
    end
    Radius_LJ_Re    = trajectInfo(2)/(2*2^(1/6));
    switch version
        case 1
            Radius_LJ_CN    = 1;
            Rbox            = trajectInfo(3);
            trajectData_3D  = zeros(Nmolecules,6,Nsamples);
        case 2
            Radius_LJ_CN    = trajectInfo(3)/(2*2^(1/6));
            Rbox            = trajectInfo(4);
            trajectData_3D  = zeros(Nmolecules,7,Nsamples);
        case 3
            Radius_LJ_CN    = trajectInfo(3)/(2*2^(1/6));
            Rbox            = trajectInfo(4);
            maxlength       = max(Nmolecules);
            trajectData_3D  = nan(maxlength,8,Nsamples);
            %%% Plot histogram/trajectory of Nmolecules
            % plot(Nmolecules); yline(mean(Nmolecules)); title(['Avg=' num2str(mean(Nmolecules)) ', shold be ' num2str(trajectInfo(1)/trajectInfo(5)) '. StDev=' num2str(std(Nmolecules))]);
            % histogram(Nmolecules); xline(mean(Nmolecules)); title(['Avg=' num2str(mean(Nmolecules)) ', shold be ' num2str(trajectInfo(1)/trajectInfo(5)) '. StDev=' num2str(std(Nmolecules))]);
    end
    
    
    if version > 2
        end_id = 0;
        for i=1:Nsamples
            start_id                = end_id + 1;
            end_id                  = start_id + Nmolecules(i) - 1;
            trajectData_3D(:,:,i)   = [trajectData_2D(start_id:end_id,:); nan(maxlength-Nmolecules(i),8)];
        end
    else
        for i=0:Nsamples-1
            start_id                    = i*Nmolecules+1;
            end_id                      = (i+1)*Nmolecules;
            trajectData_3D(:,:,i+1)     = trajectData_2D(start_id:end_id,:);
        end
    end
    
%% WRITE to dataStruct (Process 2D-IR)
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

%% WRITE to simData structure (Surface trajectories)
simData.LinSpec_FD          = LinSpec_FD;
simData.version				= version;
simData.trajectInfo			= trajectInfo;
simData.trajectData_2D		= trajectData_2D;
simData.trajectData_3D 		= trajectData_3D;
simData.Nmolecules			= Nmolecules;
simData.Nsamples			= Nsamples;
simData.Radius_LJ_Re		= Radius_LJ_Re;
simData.Radius_LJ_CN		= Radius_LJ_CN;
simData.Rbox				= Rbox;

%% Write the simData structure to dataStruct
dataStruct.simData = simData;

else
    error('Not a valid 2D-IR dataset: empty folder or corrupt data')
end