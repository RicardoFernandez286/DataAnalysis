function dataStruct = LoadDataTRUVIS(dataStruct,SpectrumToDisplay,varargin)
% Loading routine for Pump-Probe data from TRUVIS, FileFormat 2
% Updated for the new AppDesigner GUI.
%
% v4.0a / 16.04.2019 / Ricardo Fernández-Terán

dataStruct.removeduplicates = 1;
cut2ndOrder = 1;
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
    % In Lab 4, the piezo modulation is slow, so it is also included (in Lab 1 this would return -1)
    ModPiezo        = fgetl(slowmodID);
    % Close the file
fclose(slowmodID); clear slowmodID;

% Count the number of modulations for A,B,C and spectrometer:
N_modA          = sum(~isnan(str2num(ModA{2})),2);
N_modB          = sum(~isnan(str2num(ModB{2})),2);
N_modC          = sum(~isnan(str2num(ModC{2})),2);
N_modSpec       = sum(~isnan(str2num(ModSpec{2})),2);
% If any of them is zero (i.e. all NaN's), make it equal to 1 (i.e. 1 modulation = no modulation)
AllSlowMod                          = [N_modA N_modB N_modC N_modSpec];
AllSlowMod(AllSlowMod==0)           = 1;
N_modA          = AllSlowMod(1);
N_modB          = AllSlowMod(2);
N_modC          = AllSlowMod(3);
N_modSpec       = AllSlowMod(4);

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
else
    anisotropy = 'NONE';
end

% Count the number of scans
tempdir         = [rootdir filesep datafilename filesep 'temp'];
tempdirOUT      = [rootdir filesep datafilename 'temp'];
% Check if the temp dir is outside (Lab 1)
if exist(tempdir,'dir') == 0
   tempdir      = tempdirOUT;
end
tempfilelist    = dir(tempdir);
tempfilenames   = {tempfilelist.name};
% Select only "signal" files - Count the (total) number of datastates
Nscans          = round(length(tempfilenames(contains(tempfilenames,"signal",'IgnoreCase',1)))/(3.*Nslowmod));

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
cmprobe         = csvread([filename '_wavelengths.csv']);

% Read the actual data
rawsignal={}; noise={};
for m=1:Nslowmod
    for k=1:Nspectra
        ending            = ['_sp' num2str(k-1) '_sm' num2str(m-1) '_du0.csv'];
        tempsignal{m,k}   = csvread([filename '_signal' ending]);
        tempnoise{m,k}    = csvread([filename '_signal_noise' ending]);
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
end

% If Nspec = 0, then we have ONE spectrometer state (i.e. the spectrometer doesn't move)
% Otherwise, we have Nspec spectrometer states
% Recalculate the wavenumber axis according to which state is read
if N_modSpec >= 2
    WL_nm       = 10^7./cmprobe;
    NewCentreWL = SpectrometerWL(SpectrumToDisplay);
    ShiftWL     = NewCentreWL-WL_nm(floor(length(WL_nm)/2));
    cmprobe     = 10^7./(WL_nm + ShiftWL);
end

% Determine whether the electronic delay was used (ns) or the delay stage (fs)
fid         = fopen([filename '_meta.txt']); 
metadata    = textscan(fid,'%s','delimiter','\n');
metadata    = metadata{:};
fclose(fid);
ED_con_LN   = find(~cellfun(@isempty,strfind(metadata,'Electronic delay:')));
ED_con_txt  = strsplit(metadata{ED_con_LN(1)},': ');
ED_con_TF   = ~contains(ED_con_txt{2},'Not','IgnoreCase',1);
UV_con_LN   = find(~cellfun(@isempty,strfind(metadata,'UV stage:')));
UV_con_txt  = strsplit(metadata{UV_con_LN(1)},': ');
UV_con_TF   = ~contains(UV_con_txt{2},'Not','IgnoreCase',1);
Stage_LN    = find(~cellfun(@isempty,strfind(metadata,'Stage: ')));
if isempty(Stage_LN)
    if ED_con_TF && ~UV_con_TF
        % If both are connected, assume it is the EDelay
        dataStruct.timescale = 'ns';
    else
        % Otherwise, assume it is the UV stage
        delays = delays/1000;
        dataStruct.timescale = 'ps';
    end
else
    Stage_txt   = strsplit(metadata{Stage_LN(1)},'Stage: ');
    switch Stage_txt{2}
        case 'Elec delay'
            dataStruct.timescale = 'ns';
        case 'UV delay'
            delays = delays/1000;
            dataStruct.timescale = 'ps'; 
    end
end

% In TRUVIS, the 2nd order can overlap with the main band, in which case it
% should be removed (only if using a grating)
minwl = min(cmprobe);
maxwl = max(cmprobe);
if maxwl >= 2*minwl && cut2ndOrder
    NewMaxwl        = 2*minwl;
    NewMaxwl_index  = findClosestId2Val(cmprobe,NewMaxwl);
    rawsignal       = rawsignal(:,1:NewMaxwl_index);
    noise           = noise(:,1:NewMaxwl_index);
    cmprobe         = cmprobe(1:NewMaxwl_index);
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
plotranges      = [mintime maxtime minwl maxwl minabs maxabs];

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
dataStruct.cmprobe      = cmprobe;
dataStruct.rawsignal    = rawsignal;
dataStruct.noise        = noise;
dataStruct.plotranges   = plotranges;
dataStruct.Nscans       = Nscans;
dataStruct.tempdir      = tempdir;

% Calculate noise statistics
dataStruct.AvgNoise     = mean(noise(:));
dataStruct.MaxNoise     = max(noise(:));
dataStruct.SNR          = abs(round(zminmax/dataStruct.AvgNoise,3));


% Set defaults for background subtraction subunit
if dataStruct.recalcBkg == 0
    % Subtract the first negative delay from all the dataset by default - otherwise take the inputs
    dataStruct.mintimeBkg = dataStruct.delays(1);
    dataStruct.maxtimeBkg = dataStruct.delays(1);
end
    Idx = findClosestId2Val(dataStruct.delays,[dataStruct.mintimeBkg dataStruct.maxtimeBkg]);
    % Do the background subtraction and change status in handles.rawcorr
    if Idx(2)==1
        dataStruct.bkg = dataStruct.rawsignal(1,:);
    else
        dataStruct.bkg = mean(dataStruct.rawsignal(Idx(1):Idx(2),:));
    end
dataStruct.corrdata = dataStruct.rawsignal - dataStruct.bkg;
    
else
    assert(exist([filename '_delays.csv'],'file') ~= 0,'Selected directory does not contain a valid Pump-Probe (Lab 1 or Lab 4) dataset')
end

