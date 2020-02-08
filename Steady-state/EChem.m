% function EChem()
% CURRENT VERSION: v1.1a
% 04.12.2019 / Ricardo Fernandez-Teran
% Changelog:
% - Added automatic peak detection and possibility to correct errors, as well as automated Randles-Sevcik and Trumpet plots.
% - Added compatibility with Ivium data (XY XY XY XY format). Can discriminate in the case E I Q was saved instead of XY XY XY
% - Updated potential referencing and unit conversion
% - Other major improvements

%% Change the input units here
inI_units   = 'mA';     % Options: 'mA', 'microA' or 'A'
inE_units   = 'V';      % Options: 'mV' or 'V'
DataType    = 'Ivium';  % Options: 'EC-lab' or 'Ivium'

%% Sort by scan rate?
SortScanRate = 1;

%% Plotting options
WhichScan   = 2;
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
                data(:,2*j-1)   = data(:,2*j-1)*outE_factor - ref_pot_mV;
                data(:,2*j)     = data(:,2*j)*outI_factor;
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
    ScanRates_num       = ScanRates_num(sort_id(sort_SR > 0 & ~isnan(sort_SR)));
    Nfiles              = size(AllData,1);
end

%% Plot the data
fh = figure();
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

%%% Format figure and plots
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

%% Do the peak analysis
if SortScanRate ~= 1
    return
end

doAnalysis = questdlg('Do you want to perform scan-rate dependence analysis?','Scan-rate analysis','Yes','No','No');
if isempty(doAnalysis) || strcmp(doAnalysis,'No') 
    return
end
pkCurr          = zeros(Nfiles,2);
pkPot           = zeros(Nfiles,2);
E12             = zeros(Nfiles,1);
fw_E            = cell(Nfiles,1);
fw_I            = cell(Nfiles,1);
bw_E            = cell(Nfiles,1);
bw_I            = cell(Nfiles,1);
fitrange        = inputdlg('Select range for peak analysis (start,end in mV):','Select range...');
if isempty(fitrange)
    return
end
fitrange        = str2num(fitrange{1});

for i=1:Nfiles
    % First: Separate the scan in forward/backward halves
    AllE            = AllData{i,plotScan}(:,1);
    AllI            = AllData{i,plotScan}(:,2);
    InRange         = AllE>=min(fitrange) & AllE<=max(fitrange);
    AllE            = AllE(InRange);
    AllI            = AllI(InRange);
    [~,idx_max]     = max(AllE);
    [~,idx_min]     = min(AllE);
    Npoints         = length(AllE);
    idx             = 1:Npoints;
    if idx_max > idx_min % if min comes first then 1st half is reductive
        fw_pts      = [1:idx_min (idx_max+1):Npoints];
        bw_pts      = (idx_min+1):idx_max;
    else
        fw_pts      = [1:idx_max (idx_min+1):Npoints];
        bw_pts      = (idx_max+1):idx_min;
    end
    slope           = zeros(Npoints,1);
    slope(fw_pts)   = +1;
    slope(bw_pts)   = -1;
    % Second: Fit peaks on the forward/backward scans
    %%% Forward scan (reductive)
    [fw_E{i},idx_fw_E]  = sort(AllE(slope>0));
    [fw_E{i},un_fw_E]   = unique(fw_E{i});          
    fw_I{i}             = AllI(slope>0);
    fw_I{i}             = fw_I{i}(idx_fw_E(un_fw_E));
    Norm_fw             = min(fw_I{i});
    [pk_I,pk_E]         = findpeaks(fw_I{i}/Norm_fw,fw_E{i},'MinPeakHeight',0.1,'MinPeakDistance',50,'NPeaks',1);
    if ~isempty(pk_I)
        pkCurr(i,1)     = pk_I*Norm_fw;
        pkPot(i,1)      = pk_E;
    else
        pkCurr(i,1)     = 0;
        pkPot(i,1)      = 0;
    end
    %%% Backward scan (oxidative)
    [bw_E{i},idx_bw_E]  = sort(AllE(slope<0));
    [bw_E{i},un_bw_E]   = unique(bw_E{i});
    bw_I{i}             = AllI(slope<0);
    bw_I{i}             = bw_I{i}(idx_bw_E(un_bw_E));
    Norm_bw             = max(bw_I{i});
    [pk_I,pk_E]         = findpeaks(bw_I{i}/Norm_bw,bw_E{i},'MinPeakHeight',0.1,'MinPeakDistance',50,'NPeaks',1);
    if ~isempty(pk_I)
        pkCurr(i,2)     = pk_I*Norm_bw;
        pkPot(i,2)      = pk_E;
    else
        pkCurr(i,2)     = 0;
        pkPot(i,2)      = 0;
    end
    %%% E1/2
    E12(i)          = mean(pkPot(i,:));
end

%% Add peak positions to main plot
hold(ax,'on')
for i=1:Nfiles
    plot(pkPot(i,1),pkCurr(i,1),'o','Color',cmap(i,:),'HandleVisibility','off');
    plot(pkPot(i,2),pkCurr(i,2),'o','Color',cmap(i,:),'HandleVisibility','off');
end
hold(ax,'off')

%% Make Trumpet plot
fh_tr = figure('Name','Trumpet plot','NumberTitle','off');
fh_tr.Color = 'w';
ax_tr = axes('parent',fh_tr);
ax_tr.FontSize = 14;
box(ax_tr,'on');
hold(ax_tr,'on')
plot(ax_tr,ScanRates_num,pkPot(:,2)-E12,'^r','MarkerSize',8,'DisplayName','Anodic peak (oxidation)')
plot(ax_tr,ScanRates_num,pkPot(:,1)-E12,'vb','MarkerSize',8,'DisplayName','Cathodic peak (reduction)')
yline(0,'HandleVisibility','off');
legend(ax_tr,'show');
legend(ax_tr,'boxoff');
legend(ax_tr,'location','northoutside','Orientation','horizontal')
xlabel(ax_tr,'Scan rate (mV s^{-1})','FontWeight','bold')
ylabel(ax_tr,'E_{peak} - E_{1/2} (mV)','FontWeight','bold')

hold(ax_tr,'off')
ax_tr.XScale = 'log';

%% Make Randles-Sevcik plot
fh_rs = figure('Name','Randles-Sevcik plot','NumberTitle','off');
fh_rs.Color = 'w';
ax_rs = axes('parent',fh_rs);
ax_rs.FontSize = 14;
box(ax_rs,'on');
hold(ax_rs,'on')
plot(ax_rs,sqrt(ScanRates_num),abs(pkCurr(:,2)),'^r','MarkerSize',8,'DisplayName','Anodic current (oxidation)')
plot(ax_rs,sqrt(ScanRates_num),abs(pkCurr(:,1)),'vb','MarkerSize',8,'DisplayName','Cathodic current (reduction)')
legend(ax_rs,'show');
legend(ax_rs,'boxoff');
legend(ax_rs,'location','northoutside','Orientation','horizontal')
xlabel(ax_rs,'(Scan rate)^{1/2} (mV s^{-1})^{1/2}','FontWeight','bold')
ylabel(ax_rs,['Peak current (' outI_units ')'],'FontWeight','bold')
hold(ax_rs,'off')

%% Make ipc/ipa plot
fh_rv = figure('Name','Reversibility plot','NumberTitle','off');
fh_rv.Color = 'w';
ax_rv = axes('parent',fh_rv);
ax_rv.FontSize = 14;
box(ax_rv,'on');
hold(ax_rv,'on')
plot(ax_rv,ScanRates_num,abs(pkCurr(:,1)./pkCurr(:,2)),'or','MarkerSize',8)
xlabel(ax_rv,'Scan rate (mV s^{-1})','FontWeight','bold')
ylabel(ax_rv,'i_{pc}/i_{pa}','FontWeight','bold')
hold(ax_rv,'off')
ax_rv.XScale = 'log';

return;

%% Custom points (Fix when the program gets the wrong potential)
% Just need to find the right potential, it will get the current automatically
% index: (a,b)   a=Which file;  b=1 for Cathodic (negative), b=2 for Anodic (positive) peaks
a       = 3;            %#ok<*UNRCH>
b       = 1;
NewVal  = 30.08; % New X value for the min/max that we're changing

%%% Do not edit below this line
switch b
    case 1
        Xval = fw_E{a};
        Yval = fw_I{a};
    case 2
        Xval = bw_E{a};
        Yval = bw_I{a};
end
idx_new = findClosestId2Val(Xval,NewVal);

pkPot(a,b) = NewVal;
pkCurr(a,b)= Yval(idx_new);
