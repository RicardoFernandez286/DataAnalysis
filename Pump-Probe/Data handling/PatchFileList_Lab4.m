function handles = PatchFileList_Lab4(handles,sel_datatype)

% Ricardo Fernández-Terán
% v1.0 - 09.01.2017
% 
% This function will patch all the files from IR pump-probe experiments measured in Lab 4. 
% The input and output are the same Handles structure from DataAnalysis_GUI.
% 
% INPUTS:
%     handles = structure from DataAnalysis_GUI
%     sel_datatype = '1D' or '2D'
%     
% OUTPUTS: handles.$varname. Data are selected according to "sel_datatype"
%     filenames = Cell array containing the Data names of the files that contain the selected data
%     Nscans = Array containing 
%     Scanindices = Scan indices in the global list in the format [start end], for selected data

% %% DEBUG
% sel_datatype='1D';
% rootdir='D:\Ricardo Data\switchdrive\Ph.D. UZH\RESULTS\Transient absorption\Lab 4 - Pump-probe ATR\12.01.18';

%% READ from handles
rootdir = handles.rootdir;

%% Read the SampleList file
readdata = tdfread([rootdir filesep 'SampleList.txt'],'tab');

%% Parse the information from "SampleList.txt" in the chosen folder
% Convert "Type" to string array and "SampleName" to cell array
    Type=string(readdata.Type);
    SampleName=cellstr(readdata.SampleName);

% Get the scan indices for all files according to the number of scans given in SampleList
    Nscans = [1;readdata.Nscans];
    scanindices = [cumsum(Nscans(1:end-1)),cumsum(Nscans(1:end-1))+Nscans(2:end)-1];
    Nscans = Nscans(2:end);

% Select the information ONLY for the files that have the selected datatype
    SampleName = SampleName(Type == sel_datatype);
    fileindices = find(Type == sel_datatype);
    scanindices = scanindices(Type == sel_datatype,:);
    Nscans = Nscans(Type == sel_datatype);
    
%% Build the file list
    % Get a list of all files and folders in this folder.
    filelist = dir(rootdir);
    % Extract only those that are NOT directories
    filelist = filelist(not([filelist.isdir]));
    % Now extract only those that contain ".dat" in the filename (extension!)
    filelist = filelist(contains({filelist.name},'.dat'));
    % From those, extract the ones that contain "2D SS"
    filelist = filelist(contains({filelist.name},'2D SS'));
    % Create a cell array of filenames
    FileNames = {filelist.name};

%% Sort the file list according to the name format
    % NOT IMPLEMENTED!
    
    
% Not sure if the bottom sections should be here or in LoadIR4 and Load2DIR4
%% Read time zero
    delayzerofile = [rootdir filesep 'DelayZero.txt'];
    if exist(delayzerofile,'file') ~= 0
        DelayZero = csvread(delayzerofile);
    else
        DelayZero = zeros(size(scanindices,2),1);
    end
    
%% Read all the files and save all the scan data (to operate on it later)
    scandata=[];
    alldelays=[];
    tempscandata=[];
    cd(rootdir)
    for i=1:length(FileNames)
        tempscandata = dlmread(FileNames{i},'\t');
        alldelays(:,i) = tempscandata(:,1); % Time delays in the first column
        scandata(:,:,i) = tempscandata(:,2:end); % DeltaAbs (mOD) in second column of data
    end

%% Verify if "Wavenumbers.dat" exists and read it, otherwise plot with pixel numbers
    if exist('Wavenumbers.dat','file') ~= 0
        cmprobe = dlmread('Wavenumbers.dat','\t');
        WLcalibrated = 1;
    else
        cmprobe = 1:1:size(scandata,2);
        WLcalibrated = 0;
    end
    
%% WRITE to handles
    % Filenames
    handles.SampleName = SampleName;
    handles.filenames = FileNames;
    handles.fileindices = fileindices;
    handles.scanindices = scanindices;
    handles.Nscans = Nscans;
    % Data
    handles.alldelays = alldelays;  % Double (Delays x scans)
    handles.cmprobe = cmprobe;      % Double (pixels x 1)
    handles.scandata = scandata;    % Double (Delays x pixels x scans)
    handles.DelayZero = DelayZero;  % Double (1 x Nsamples)