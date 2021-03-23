function ALLdata = batchCLS(rootfolder,Nskip,fileList,varargin)
%% Define startup variables
% rootfolder  = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\9) 2D IR distance - na\Latest Simulations\Dimer_distance1';
% fitType     = 'Test'; % 'Ricardo' or 'Andrea'

%% Build list of (sub)folders, which contain the data
folderslist = dir(rootfolder);
foldernames = {folderslist.name}';
foldernames = foldernames([folderslist.isdir]);
foldernames = foldernames(3:end);
Ndatafiles  = length(foldernames);

% Read first options
if ~isa(Nskip,'double')
    Nskip = 1;
end
disp([datestr(now,-1) ': Fitting every ' num2str(Nskip) ' points in t2.']);

if isa(fileList,'string') && ~isempty(fileList)
    disp([datestr(now,-1) ': Reading directory list...']);
    fid = fopen(fileList);
    j=0;
    while ~feof(fid)
        j=j+1;
        imp_datapath{j} = fgetl(fid);
    end
    imp_datapath = imp_datapath';
    Ndatafiles = length(imp_datapath);
    fclose(fid);
    disp([datestr(now,-1) ': Found ' num2str(Ndatafiles) ' datasets in list.']);
end

% Quit if nothing to do
if Ndatafiles == 0
    disp([datestr(now,-1) ': Nothing to do here! Bye :)']);
%     exit
    return
end

%% Loop through the folders: Load, process and analyse all the spectra
for i=1:Ndatafiles
    if length(varargin) > 1
%         [folderpath,foldernames{i}] = fileparts(imp_datapath{i});
        imp_parts   = strsplit(imp_datapath{i},filesep);
        folderpath  = [];
        for s=1:length(imp_parts)-1
            folderpath = [folderpath filesep imp_parts{s}];
        end
        foldernames{i} = imp_parts{end};
        rootfolder  = ['/home/group' filesep folderpath];
    end

    %%% Determine the data type
    if exist([rootfolder filesep foldernames{i} filesep 'spec_2D.dat'],'file') ~= 0
        dataType                = 'Simulation';
    elseif exist([rootfolder filesep foldernames{i} filesep foldernames{i} '_bins.csv'],'file') ~= 0
        dataType                = 'Experiment';
    else
        disp([datestr(now,-1) ': ' 'Folder ' foldernames{i} ' does not contain a valid 2D-IR dataset']);
        continue % Not a valid dataset, try with the the next one
    end
    
    %%% Load the data
    dataStruct.rootdir      = rootfolder;
    dataStruct.datafilename = foldernames{i};
    load([fileparts(mfilename('fullpath')) filesep 'defaultGUI.mat']);
    
    switch dataType
        case 'Simulation'
            disp([datestr(now,-1) ': ' 'Loading dataset ' foldernames{i} ' [Simulation]']);
            dataStruct = load2DIRsimu(dataStruct,'NoWaitBar');
        case 'Experiment'
            disp([datestr(now,-1) ': ' 'Loading dataset ' foldernames{i} ' [Experiment]']);
            dataStruct.Transient2D = 0;
            dataStruct = load2DIRlab1(dataStruct,'NoWaitBar');
    end
    
%% Define important parameters
t2_fitrange         = [0.25 30];

% writepath           = '/home/ricfer/FitResults/Solvation/IrHCOP3';

writepath = [rootfolder filesep 'CLS'];

intensity_threshold = 50/100;
Interp_method       = 'InterpFT';

%%%%% FOR IrHCOP3
% pump_ranges         = [1900 1950; 2050 2100]; % Each row is a range
% probe_ranges        = [1900 1950; 2050 2100]; % Each row is a range

%%%%% FOR VC-H2 low freq
% pump_ranges         = [1970 2020; 2070 2120]; % Each row is a range
% probe_ranges        = [1970 2010; 2070 2120]; % Each row is a range

%%%%% FOR VC-H2 high freq
pump_ranges         = [2070 2120; 2170 2260]; % Each row is a range
probe_ranges        = [2070 2120; 2180 2260]; % Each row is a range


%% Write settings to structure
settings.intensity_threshold    = intensity_threshold;
settings.Interp_method          = Interp_method;
settings.pump_ranges            = pump_ranges;
settings.probe_ranges           = probe_ranges;
settings.t2_fitrange            = t2_fitrange;

%% Load the data and do the analysis  
% if exist([writepath filesep foldernames{i} '_CLS.csv'],'file') ~= 0
%     disp([datestr(now,-1) ': ' 'Dataset ' foldernames{i} ' has been analysed before! Delete/rename the old analysis file and try again.']);
%     continue
% end

    %%% Process the data (experiment only)
    switch dataType
        case 'Experiment'
            disp([datestr(now,-1) ': ' 'Processing ' foldernames{i} ' ...']);
            app.I2D_AutocalibrateprobeaxisCheckBox.Value = 0;
            app.I2D_ZeropadFactor.Value = 4;
            app.I2D_ApodisationFunction.Value = '2'; % Cos^n apodisation  n=0 for box.
            dataStruct = process2DIR(app,dataStruct,0,'NoWaitBar');
    end

    %%% Do the spectral diffusion analysis
    try
        disp([datestr(now,-1) ': ' 'Starting CLS analysis of ' foldernames{i} ' ...']);
        dataStruct      = SpectralDiffusion_analysis(app,dataStruct,[],'Batch',settings);
        dataOUT         = [dataStruct.AnalysedDelays dataStruct.SpecDif_ind_all];
        outFileName     = [writepath filesep foldernames{i} '_CLS.csv'];
        if exist(writepath,'dir') == 0; mkdir(writepath); end
        writematrix(dataOUT,outFileName);
        disp([datestr(now,-1) ': Analysis of dataset  ' foldernames{i} '  is  [DONE]']);
        ALLdata(:,:,i)  = dataOUT;
    catch
        disp([datestr(now,-1) ': ' 'Error analysing ' foldernames{i} ' ...']);
        continue
    end
end

disp([datestr(now,-1) ': ' 'Everything done!']);


plotMHsolv_CLS;
% exit