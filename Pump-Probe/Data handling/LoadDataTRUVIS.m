function dataStruct = LoadDataTRUVIS(dataStruct,SpectrumToDisplay,varargin)
% Loading routine for Pump-Probe data from TRUVIS, FileFormat 2
% Updated for the new AppDesigner GUI.
% Updated for compatibility with n-detector arrangement
%
% v5.0b / 25.05.2022 / Ricardo Fernández-Terán

dataStruct.removeduplicates = 1;
cut2ndOrder = 0;
%% READ from dataStruct
rootdir         = dataStruct.rootdir;
datafilename    = dataStruct.datafilename;

%% Prepare the filenames for each file to read
filename        = [rootdir filesep datafilename filesep datafilename];

%% Check that all files exist and read them

% Read the files
if exist([filename '_delays.csv'],'file') ~= 0 % If all files exist
% COMMENT:
%   SP = Fast modulation spectra
%   SM = Slow modulation (as defined in the SlowMod cluster)

% Consider the datastates
DatastateMatrix = csvread([filename '_spectraAndDatastates.csv']);
Ndatastates     = length(unique(DatastateMatrix(:,2)));
Nspectra        = length(unique(DatastateMatrix(:,1)));

% Read all the lines of the SlowModulation file
slowmodID       = fopen([filename '_slowModulation.csv']);
    % The SlowModulation file contains the No. of slowmods in the first line
    Nslowmod        = str2double(fgetl(slowmodID));
    % The rest of the lines contain the action matrices of all different modulations
    ModA            = strsplit(fgetl(slowmodID),':');
    ModB            = strsplit(fgetl(slowmodID),':'); % IR pol
    ModC            = strsplit(fgetl(slowmodID),':'); % UV pol
    ModSpec         = strsplit(fgetl(slowmodID),':');
    
    %% LAB 4
    % % In Lab 4, the piezo modulation is slow, so it is also included (in Lab 1 this would return -1)
    % ModPiezo        = fgetl(slowmodID);
    ModPiezo = -1;
    %% UniGE
    % At UniGE, the third modulation are the waveplates
    ModWP = strsplit(fgetl(slowmodID),':');

    % Close the file
fclose(slowmodID); clear slowmodID;

% Count the number of modulations for A,B,C and spectrometer:
N_modA          = sum(~isnan(str2num(ModA{2})),2);
N_modB          = sum(~isnan(str2num(ModB{2})),2);
N_modC          = sum(~isnan(str2num(ModC{2})),2);
N_modSpec       = sum(~isnan(str2num(ModSpec{2})),2);
N_modWP         = sum(~isnan(str2num(ModWP{2})),2);
% If any of them is zero (i.e. all NaN's), make it equal to 1 (i.e. 1 modulation = no modulation)
AllSlowMod                          = [N_modA N_modB N_modC N_modSpec N_modWP];
AllSlowMod(AllSlowMod==0)           = 1;
N_modA          = AllSlowMod(1);
N_modB          = AllSlowMod(2);
N_modC          = AllSlowMod(3);
N_modSpec       = AllSlowMod(4);
N_modWP         = AllSlowMod(5);
% Get the spectrometer wavelengths
SpectrometerWL      = str2num(ModSpec{2});

% Count the number of 4-point modulation states (Lab 4 only)
if ModPiezo ~= -1 % if not at the end of the file
    ModPiezo    = strsplit(ModPiezo,':');
    N_4PM       = length(str2num(ModPiezo{2}))./(N_modA*N_modB*N_modC*N_modSpec);
else
    N_4PM = 1;
end

% Consider anisotropy or magic angle signal calculation
if ~isempty(varargin)
    anisotropy = varargin{1};
    recalcscans = varargin{2};
else
    anisotropy = 'NONE';
    recalcscans = 0;
end

% Count the number of scans
tempdir         = [rootdir filesep datafilename filesep 'temp'];
tempdirOUT      = [rootdir filesep datafilename 'temp'];
Nscans          = NaN;
% Check if the temp dir is outside (Lab 1)
if exist(tempdir,'dir') == 0
    if exist(tempdirOUT,'dir') ~= 0
        tempdir = tempdirOUT;
        Nscans  = 0;
    end
else
    Nscans = 0;
end

if ~isnan(Nscans)
    tempfilelist    = dir(tempdir);
    tempfilenames   = {tempfilelist.name};
    % Select only "signal" files - Count the (total) number of datastates
    Nscans          = round(length(tempfilenames(contains(tempfilenames,"signal",'IgnoreCase',1)))/(3.*Nslowmod));
end

% Whatever is not a 4-point modulation is a real "modulation"
mod             = Nslowmod/N_4PM;
SM_selString    = cell(mod,1);
for i=1:mod
    SM_selString{i} = num2str(i);
end
dataStruct.SlowMod_selector.String     = SM_selString;

%% Define the measurement scheme
switch Ndatastates
    case 1
        chopper = 'OFF';
    case 2
        chopper = 'ON';
end

switch N_4PM
    case 1
        piezomod = 'OFF';
    case 4
        piezomod = '4PM';
end

delays          = csvread([filename '_delays.csv']);
% If there is a calibrated probe file, use it. Otherwise, use the normal one
if exist([rootdir filesep datafilename filesep 'CalibratedProbe.csv'],'file') == 2
    cmprobe     = csvread([rootdir filesep datafilename filesep 'CalibratedProbe.csv']);
elseif exist([rootdir filesep 'CalibratedProbe.csv'],'file') == 2
    cmprobe     = csvread([rootdir filesep 'CalibratedProbe.csv']);
else
    cmprobe     = csvread([filename '_wavelengths.csv']);
end

% Read the actual data
rawsignal   = {};
noise       = {};
tempsignal  = cell(Nslowmod,Nspectra);
tempnoise   = cell(Nslowmod,Nspectra);

for m=1:Nslowmod
    for k=1:Nspectra
        ending            = ['_sp' num2str(k-1) '_sm' num2str(m-1) '_du0.csv'];
        tempsignal{m,k}   = readmatrix([filename '_signal' ending],'FileType','text');
        tempnoise{m,k}    = readmatrix([filename '_signal_noise' ending],'FileType','text');
    end
end

switch anisotropy
    case 'NONE'
        switch piezomod
            case 'OFF' % The signal is the signal, just read it from the cell array
                rawsignal         = tempsignal{SpectrumToDisplay,1};
                noise             = tempnoise{SpectrumToDisplay,1};
            case '4PM' % The signal comes from a four point modulation, average the four points
                rawsignal         = sum(cat(3,tempsignal{:}),3)/4;
                noise             = sum(cat(3,tempnoise{:}),3)/4;
        end
    case 'UV anisotropy'
        par               = tempsignal{2,1};
        perp              = tempsignal{1,1};
        par_n             = tempnoise{2,1};
        perp_n            = tempnoise{1,1};
        
        rawsignal         = (par - perp)./(par + 2.*perp);
        noise             = abs(((par+2.*perp)-(par-perp))./(par+2.*perp).^2).*par_n + abs((-(par+2.*perp)-2.*(par-perp))./(par+2.*perp).^2).*perp_n;
    case 'IR anisotropy'
        par               = tempsignal{1,1};
        perp              = tempsignal{2,1};
        par_n             = tempnoise{1,1};
        perp_n            = tempnoise{2,1};
        
        rawsignal         = (par - perp)./(par + 2.*perp);
        noise             = abs(((par+2.*perp)-(par-perp))./(par+2.*perp).^2).*par_n + abs((-(par+2.*perp)-2.*(par-perp))./(par+2.*perp).^2).*perp_n;
    case 'UV magic angle'
        par               = tempsignal{2,1};
        perp              = tempsignal{1,1};
        par_n             = tempnoise{2,1};
        perp_n            = tempnoise{1,1};
        
        rawsignal         = (par + 2.*perp)./3;
        noise             = (par_n + 2.*perp_n)./3;
    case 'IR magic angle'
        par               = tempsignal{1,1};
        perp              = tempsignal{2,1};
        par_n             = tempnoise{1,1};
        perp_n            = tempnoise{2,1};
        
        rawsignal         = (par + 2.*perp)./3;
        noise             = (par_n + 2.*perp_n)./3;
    case 'UniGE Magic Angle (calc.)'
        par               = tempsignal{1,1};
        perp              = tempsignal{3,1};
        par_n             = tempnoise{1,1};
        perp_n            = tempnoise{3,1};
        
        par = par - mean(par(1:3,:));
        perp = perp - mean(perp(1:3,:));

        rawsignal         = (par + 2.*perp)./3;
        noise             = (par_n + 2.*perp_n)./3;
    case 'UniGE Magic Angle (exp.)'
        MA                = tempsignal{2,1};
        
        rawsignal         = MA - mean(MA(1:3,:));
        noise             = tempnoise{2,1};


    case 'UniGE Anisotropy'
        par               = tempsignal{1,1};
        perp              = tempsignal{3,1};
        par_n             = tempnoise{1,1};
        perp_n            = tempnoise{3,1};
        
        par = par - mean(par(1:3,:));
        perp = perp - mean(perp(1:3,:));

        rawsignal         = (par - perp)./(par + 2.*perp);
        noise             = abs(((par+2.*perp)-(par-perp))./(par+2.*perp).^2).*par_n + abs((-(par+2.*perp)-2.*(par-perp))./(par+2.*perp).^2).*perp_n;
    case 'UniGE MA check'
        par               = tempsignal{1,1};
        perp              = tempsignal{3,1};
        par_n             = tempnoise{1,1};
        perp_n            = tempnoise{3,1};
        
        par = par - mean(par(1:3,:));
        perp = perp - mean(perp(1:3,:));

        MA                = tempsignal{2,1};
        MA_n              = tempnoise{2,1};
      
        MA = MA - mean(MA(1:3,:));

        rawsignal         = MA - (par + 2.*perp)./3;
        noise             = MA_n + (par_n + 2.*perp_n)./3;
    case 'UniGE Pol. diff. (Par-Perp)'
        rawsignal         = tempsignal{1,1}-tempsignal{3,1};
        noise             = tempnoise{1,1}+tempnoise{3,1};
end

if recalcscans == 1
    rawsignal       = dataStruct.rawsignal;
    noise           = dataStruct.noise;
end

% If there is a calibrated probe file, use it. Otherwise, use the normal one
if exist([rootdir filesep datafilename filesep 'CalibratedProbe.csv'],'file') == 2
    cmprobe     = csvread([rootdir filesep datafilename filesep 'CalibratedProbe.csv']);
elseif exist([rootdir filesep 'CalibratedProbe.csv'],'file') == 2
    cmprobe     = csvread([rootdir filesep 'CalibratedProbe.csv']);
else
    cmprobe     = csvread([filename '_wavelengths.csv']);
end

% In TRUVIS, the 2nd order can overlap with the main band, in which case it
% should be removed (only if using a grating)
minwl = min(cmprobe);
maxwl = max(cmprobe);

if unique(diff(cmprobe))==1 
    cut2ndOrder = 0;
end

if maxwl >= 2*minwl && cut2ndOrder
    NewMaxwl        = 2*minwl;
    NewMaxwl_index  = findClosestId2Val(cmprobe,NewMaxwl);
    rawsignal       = rawsignal(:,1:NewMaxwl_index);
    noise           = noise(:,1:NewMaxwl_index);
    cmprobe         = cmprobe(1:NewMaxwl_index);
end

% This needs a rework, but can't think of an easy way
% if the negative time delay is <10 ps = <10000 fs then we know we are in ps
if log10(abs(delays(1))) >= 3
    delays = delays/1000;
    dataStruct.timescale = 'ps';
    % Timescale is in picoseconds! (nicer to understand)
else
    dataStruct.timescale = 'ns';
end

% Read the plot ranges
mintime         = -0.5;
maxtime         = max(delays);
flatbl          = findClosestId2Val(delays,mintime);
minabs          = min(rawsignal(flatbl:end,:),[],'all');
maxabs          = max(rawsignal(flatbl:end,:),[],'all');
zminmax         = round(max([abs(minabs) abs(maxabs)]),3);
minwl           = min(cmprobe);
maxwl           = max(cmprobe);
Ncontours       = 40;
plotranges      = [mintime maxtime minwl maxwl minabs maxabs Ncontours];

% Remove duplicate delays and average the data for repeated delays (new universal method)
[delays,~,idx]  = unique(delays,'stable');
new_rawsignal   = zeros(length(delays),size(rawsignal,2));
new_noise       = zeros(length(delays),size(rawsignal,2));

for i=1:size(rawsignal,2)
    new_rawsignal(:,i)       = accumarray(idx,rawsignal(:,i),[],@mean); 
    new_noise(:,i)           = accumarray(idx,noise(:,i),[],@mean); 
end

% Sort the delays from negative to positive
[delays,sortID] = sort(delays);
rawsignal       = new_rawsignal(sortID,:);
noise           = new_noise(sortID,:);

%% Write to DataStructure
dataStruct.delays       = delays;
dataStruct.cmprobe      = {cmprobe};
dataStruct.rawsignal    = {rawsignal};
dataStruct.noise        = {noise};
dataStruct.plotranges   = {plotranges};
dataStruct.Nscans       = Nscans;
dataStruct.tempdir      = tempdir;

% Calculate noise statistics
dataStruct.AvgNoise     = mean(noise(:));
dataStruct.MaxNoise     = max(noise(:));
dataStruct.SNR          = abs(round(zminmax/dataStruct.AvgNoise,3));

switch anisotropy
    case 'NONE'
        % Set defaults for background subtraction subunit
        if dataStruct.recalcBkg == 0
            % Subtract the first negative delay from all the dataset by default - otherwise take the inputs
            dataStruct.mintimeBkg = dataStruct.delays(1);
            dataStruct.maxtimeBkg = dataStruct.delays(1);
        end
            Idx = findClosestId2Val(dataStruct.delays,[dataStruct.mintimeBkg dataStruct.maxtimeBkg]);
            % Do the background subtraction and change status in handles.rawcorr
            if Idx(2)==1
                dataStruct.bkg = {dataStruct.rawsignal{1}(1,:)};
            else
                dataStruct.bkg = {mean(dataStruct.rawsignal{1}(Idx(1):Idx(2),:))};
            end
        dataStruct.corrdata = {dataStruct.rawsignal{1} - dataStruct.bkg{1}};
    otherwise
        %do nothing
        dataStruct.corrdata = dataStruct.rawsignal;
end
    
else
    assert(exist([filename '_delays.csv'],'file') ~= 0,'Selected directory does not contain a valid Pump-Probe (Lab 1 or Lab 4) dataset')
end

