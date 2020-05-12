function batchPRDF(rootfolder,varargin)
%% Build list of (sub)folders, which contain the data
folderslist  = dir(rootfolder);
foldernames = {folderslist.name}';
foldernames = foldernames([folderslist.isdir]);
foldernames = flipud(foldernames(3:end));
Ndatafiles  = length(foldernames);

writepath   = [rootfolder filesep 'pRDF'];

if exist(writepath,'dir') == 0
    mkdir(writepath);
end

if exist([rootfolder filesep '1D'],'dir') == 0
    mkdir([rootfolder filesep '1D']);
end

% Read more options from varargin:
%   varargin{1} = No. of t2 delays to skip
%   varargin{2} = Path to a text file containing a list of the datasets to fit
if ~isempty(varargin)   
    disp([datestr(now,-1) ': Reading directory list...']);
    fid = fopen(varargin{2});
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

if Ndatafiles == 0
    disp([datestr(now,-1) ': Nothing to do here! Bye :)']);
    exit
end

%% Start parallel pool (if none exists)
% if isempty(gcp('nocreate'))
%    parpool;
% end

%% Loop through the folders: Load, process and fit all the spectra
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
            disp([datestr(now,-1) ': ' 'Skipping dataset ' foldernames{i} ' [Experiment]']);
            continue
    end
         
    if exist([writepath filesep foldernames{i} '_pRDF.dat'],'file') ~= 0
        disp([datestr(now,-1) ': ' 'Dataset ' foldernames{i} ' has been analysed before! Delete/rename the old file before trying to analyse again.']);
        continue
    end

    %%% Calculate and save pRDF
    plotRDF(dataStruct.simData,writepath,foldernames{i},'On')
    
    %%% Save 1D spectra as well
    writematrix(dataStruct.simData.LinSpec_FD,[rootfolder filesep '1D' filesep foldernames{i} '_1D.dat']);
end

disp([datestr(now,-1) ': ' 'Everything done!']);
% exit