function EChem()
% CURRENT VERSION: v1.1a
% 26.11.2019 / Ricardo Fernandez-Teran
% Changelog:
% - Added compatibility with Ivium data (XY XY XY XY format). Can discriminate in the case E I Q was saved instead of XY XY XY
% - Updated potential referencing and unit conversion
% - Other major improvements

%% Change the input units here
inI_units   = 'mA';     % Options: 'mA', 'microA' or 'A'
inE_units   = 'V';      % Options: 'mV' or 'V'
DataType    = 'Ivium'; % Options: 'EC-lab' or 'Ivium'

%% Sort by scan rate?
SortScanRate = 0;

%% Plotting options
WhichScan   = 1;
LineWidth   = 1.5;

%% Get the directory and a list of .txt files to load and plot
rootdir     = uigetdir;

if rootdir == 0
    return
end

filelist    = dir(rootdir);
filenames   = {filelist.name}';
filenames   = filenames(contains(filenames,'.txt','IgnoreCase',true));
[filenames_idx,filesSelected] = listdlg('ListString',filenames,'SelectionMode','multiple','ListSize',[400,300],'PromptString','Select files to plot and load:');

if filesSelected == 0
    return
end

ref_type    = questdlg('Which internal reference did you use?','Select internal reference','Ferrocene','Fc(Cp*)2','None/Other','Ferrocene');

switch ref_type
    case 'Ferrocene'
        ref_pot         = inputdlg('Ferrocene half-potential (mV):','Define potential scale',[1,50],{'0'});
        if(isempty(ref_pot))
            return
        end
        ref_pot_mV      = str2num(ref_pot{1});
        ref_pot_name    = 'Fc/Fc^+';
    case 'Fc(Cp*)2'
        ref_pot         = inputdlg({'Reference couple half-potential (mV):';'Potential of Ferrocene in your solvent (mV vs. Fc(Cp*)2):'},'Define potential scale',[1,50],{'0';'MeCN=505; THF=427'});
        if(isempty(ref_pot))
            return
        end
        ref_pot_mV      = str2num(ref_pot{1})+str2double(ref_pot{2});
        ref_pot_name    = 'Fc/Fc^+';
    case 'None/Other'
        ref_pot         = inputdlg({'Reference couple half-potential (mV):';'Reference couple name:'},'Define potential scale',[1,50],{'0';'Fc/Fc^+'});
        if(isempty(ref_pot))
            return
        end
        ref_pot_mV      = str2num(ref_pot{1});
        ref_pot_name    = ref_pot{2};
end


%% Define plotting units and options

% Change units
switch inI_units
    case 'mA'
        outI_factor = 1;
    case 'microA'
        outI_factor = 10^-3;
    case 'A'
        outI_factor = 10^6;
end

switch inE_units
    case 'mV'
        outE_factor = 1;
    case 'V'
        outE_factor = 1000;
end

outI_units  = 'mA';   % set the legend of the current axis

if ~isempty(ref_pot_name)
    outE_units  = ['mV vs. ' ref_pot_name];     % set the legend of the potential axis (reference comes later)
else
    outE_units  = 'mV' ;     % set the legend of the potential axis (reference comes later)
end

%% Load the files and separate into cell array CV data vs scan number
% Data structure: Cell array(File index, ScanNo) / [E x I]
Nfiles      = length(filenames_idx);
Nscans      = zeros(Nfiles,1);
AllData     = cell(Nfiles,1);
caption     = cell(Nfiles,1);

switch DataType
    case 'Ivium'
        for i=1:Nfiles
            % Read from file
            data        = readmatrix([rootdir filesep filenames{filenames_idx(i)}],'NumHeaderLines',1);
            % Remove NaN columns
            data        = data(:,any(~isnan(data), 1));
            % Sort into {[E x I]} vs scan number
            Nscans(i)   = floor(size(data,2)/2); 
            
            % Sort into {[E x I]} vs scan number
            for j=1:Nscans(i)
                data(:,1)   = data(:,1)*outE_factor - ref_pot_mV;
                data(:,2)   = data(:,2)*outI_factor;
                AllData{i,j}    = data(:,(2*j-1):2*j);
            end
            % Get a caption string
            [~,caption{i},~]    = fileparts(filenames{filenames_idx(i)});
        end
    case 'EC-lab'
        for i=1:Nfiles
            % Read from file
            data        = readmatrix([rootdir filesep filenames{filenames_idx(i)}],'NumHeaderLines',1);
            data(:,1)   = data(:,1)*outE_factor - ref_pot_mV;
            data(:,2)   = data(:,2)*outI_factor;
            % Sort into {[E x I]} vs scan number
            Nscans(i)   = max(data(:,3)); 
            for j=1:Nscans
                AllData{i,j}    = data(data(:,3)==j,1:2);
            end
            % Get a caption string
            [~,caption{i},~]    = fileparts(filenames{filenames_idx(i)});
        end
end

%% Sort scan rate dependence
if SortScanRate == 1
    ScanRates_st        = inputdlg(caption,'Input scan rate (mV/s), -1 to ignore');
    
    if isempty(ScanRates_st)
        return
    end
    Nfiles              = length(ScanRates_st);
    ScanRates_num       = zeros(Nfiles,1);
    for i=1:Nfiles
        ScanRates_num(i)= str2double(ScanRates_st{i});
        if ScanRates_num(i) > 1000
            caption{i}      = [num2str(ScanRates_num(i)/1000) ' V/s'];
        else
            caption{i}      = [num2str(ScanRates_num(i)) ' mV/s'];
        end
    end
    
    [sort_SR,sort_id]   = sort(ScanRates_num,'ascend');
    AllData             = AllData(sort_id(sort_SR > 0 & ~isnan(sort_SR)),:);
    caption             = caption(sort_id(sort_SR > 0 & ~isnan(sort_SR)));
    Nscans              = Nscans(sort_id(sort_SR > 0 & ~isnan(sort_SR)));           
    Nfiles              = size(AllData,1);
end

%% Plot
fh = figure(1);
clf(fh);
ax = axes('parent',fh);

hold(ax,'on')
if SortScanRate == 1
    cmap = colormap(ax,othercolor('Mrainbow',Nfiles));
    for i=1:Nfiles
        plotScan = min(Nscans(i),WhichScan);
        plot(AllData{i,plotScan}(:,1),AllData{i,plotScan}(:,2)*outI_factor,'DisplayName',caption{i},'Color',cmap(i,:),'LineWidth',LineWidth)
    end
else
    for i=1:Nfiles
        plotScan = min(Nscans(i),WhichScan);
        plot(AllData{i,plotScan}(:,1),AllData{i,plotScan}(:,2)*outI_factor,'DisplayName',caption{i},'LineWidth',LineWidth)
    end
end
hold(ax,'off');

%% Format figure and plots
% Figure format
fh.Color    = 'w';
fh.Position = [300,600,1160,420];

% Axis format
legend(ax,'show');
legend(ax,'location','eastoutside');
legend(ax,'interpreter','none');
legend(ax,'boxoff');
ax.Position = [0.07 0.1475 0.60 0.8];

box(ax,'on');
ax.FontSize = 14;
xlabel(ax,['Potential (' outE_units ')'],'FontWeight','bold');
ylabel(ax,['Current (' outI_units ')'],'FontWeight','bold');
axis(ax,'tight');

% Add XY lines
xl = xline(ax,0,'HandleVisibility','off'); xl.Color = 0.75*[1 1 1];
yl = yline(ax,0,'HandleVisibility','off'); yl.Color = 0.75*[1 1 1];


