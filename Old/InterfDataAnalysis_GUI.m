function varargout = InterfDataAnalysis_GUI(varargin)
%INTERFDATAANALYSIS_GUI MATLAB code file for InterfDataAnalysis_GUI.fig
%      INTERFDATAANALYSIS_GUI, by itself, creates a new INTERFDATAANALYSIS_GUI or raises the existing
%      singleton*.
%
%      H = INTERFDATAANALYSIS_GUI returns the handle to a new INTERFDATAANALYSIS_GUI or the handle to
%      the existing singleton*.
%
%      INTERFDATAANALYSIS_GUI('Property','Value',...) creates a new INTERFDATAANALYSIS_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to InterfDataAnalysis_GUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      INTERFDATAANALYSIS_GUI('CALLBACK') and INTERFDATAANALYSIS_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in INTERFDATAANALYSIS_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help InterfDataAnalysis_GUI

% Last Modified by GUIDE v2.5 18-Jul-2018 19:45:33

% Ricardo Fern�ndez-Ter�n, v4.8c - 29.03.2019

% ----CHANGELOG:
% * Fixed bug when loading raw Transient 2D data. Now it should be a bit more general.
% * Data list will now omit and not show the "temp" folders
% * Implemented spectral subtraction module
% * Implemented spectral diffusion module
% * Added compatibility with transient 2D IR datasets
% * Updated plotting routine: now the contour levels are defined based on min/max/# levels
% * Updated GUI features, including progress bar.
% * Implemented plotting of t2 kinetics and spectral slices (integrated still pending).
% * Fixed phasing routine - now the phase is flat around the spectral maximum (as it should!)
% * Fixed loading routine: now it is possible to read both SIGNAL and RAW (intensity) data files
% * Added "Go to today" button, to go to the current day's folder.

% ----PENDING:
% * Movie, actions of most of the buttons :)
% * Plot Z range control in main screen (!!)
% * Anisotropy calculations from two datasets


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @InterfDataAnalysis_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @InterfDataAnalysis_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before InterfDataAnalysis_GUI is made visible.
function InterfDataAnalysis_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for InterfDataAnalysis_GUI
handles.output = hObject;

% Set default directory:
if exist('C:\GUIoptions.txt','file')==2
    handles.defaultdir = readParam('default2DIRdir','C:\GUIoptions.txt');
else
    handles.defaultdir = [];
end

% Update version text string
handles.VersionText.String = "v4.8c - 29.03.2019";

% Disable annoying warnings
warning('off','MATLAB:Axes:NegativeLimitsInLogAxis');
warning('off','MATLAB:legend:IgnoringExtraEntries');

% UIWAIT makes InterfDataAnalysis_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = InterfDataAnalysis_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in BrowseRootDirButton.
function BrowseRootDirButton_Callback(hObject, eventdata, handles)
rootdir = uigetdir(handles.defaultdir); 
handles.rootdir = rootdir;
if rootdir ~= 0 % the user selected a directory
    % Update text label
    set(handles.CurrDir, 'String', rootdir);
    if handles.DataTypeMenu.Value ~= 4
        % Put list of subfolders into subFolders
            % Get a list of all files and folders in this folder.
            files = dir(rootdir);
            % Extract only those that are directories.
            subFolders = files([files.isdir]);
            % Don't show the temp folders in the list
            subFolders = subFolders(~contains({subFolders.name},"temp",'IgnoreCase',1));
        % Go to the selected folder
            cd (rootdir);
        % Sort the names and create variables for the list
            dir_struct = subFolders;
            [sorted_names,sorted_index] = sortrows({dir_struct.name}');
            handles.file_names = sorted_names;
            handles.is_dir = [dir_struct.isdir];
            handles.sorted_index = sorted_index;
            guidata(hObject,handles)
        % Update DatafoldersList, removing ".." and "." from the visible list
            handles.file_names = handles.file_names(3:end);
            set(handles.DatafoldersList,'String',handles.file_names,'Value',1)
            guidata(hObject,handles)
    else
        filelist            = dir(rootdir);
        filenames           = {filelist.name}';
        filenames           = filenames(~[filelist.isdir]);
        filenames_delays    = filenames(~contains(filenames,'Pinhole','IgnoreCase',true));
        filenames_delays    = filenames_delays(~contains(filenames_delays,'Phases','IgnoreCase',true));
               
        % Read the formatted strings of each filename
        Ndelays         = length(filenames_delays);
        time            = zeros(1,Ndelays);
        for i=1:Ndelays
            parts       = strsplit(filenames_delays{i},'_');
            time_str    = strsplit(parts{end-3},'fs');
            time(i)     = str2double(time_str{1})./1000;
        end
        SampleName      = [parts{2:end-4}];
        
        % WRITE TO HANDLES
        handles.t2delays    = time';
        handles.Ndelays     = Ndelays;
        handles.filenames   = filenames_delays;
        % Write to the front panel
        handles.Population_delay.String = num2str(sort(handles.t2delays));
        set(handles.DatafoldersList,'String',SampleName,'Value',1)
        handles.WaitBar         = waitbar(0,'Loading data...');
        for j=0:Ndelays
            handles             = load2DIRlab3(handles,j);
            guidata(hObject,handles)
        end
        % Delete progress bar
        delete(handles.WaitBar);
        % Update the edit controls
        handles = UpdateEdits_2DIR(handles,1);
        % Plot the data in the main window (for preview)
        plot_2DIR(handles,handles.MainAxes);   
        set(handles.MainAxes,'Visible','On')
        % Generate the secondary plot (show phasing data by default)
%         SecondaryPlot(handles,handles.SecondaryAxes,'td');
        set(handles.SecondaryAxes,'Visible','Off')
        % Enable disabled controls
        EnableControls_2DIR(handles)
        % Update handles
        guidata(hObject,handles)
    end
else 
    handles.CurrDir.String = 'Please select a directory!  :(';
end

% --- Executes on selection change in DatafoldersList.
function DatafoldersList_Callback(hObject, eventdata, handles)
if handles.DataTypeMenu.Value ~= 4 % Only do something if not Lab 3
    
% Determine whether the user double-clicked the folder or not
switch get(gcf,'SelectionType')
case 'normal' % Single click
    index_selected  = handles.DatafoldersList.Value;
    file_list       = handles.DatafoldersList.String;
    datafilename    = cellstr(file_list{index_selected}); % Item selected in list box
    % Write datafilename to handles(datafilename) and update
    handles.datafilename = string(datafilename);
    guidata(hObject,handles)
    % If directory
    if  handles.is_dir(handles.sorted_index(index_selected))
        % Load the data (replace handles) and update the main handles
        error=0;
        try
           % Check the number of t2 delays
           OldNt2delays = length(handles.Population_delay.String);
           %%% Load the data
           switch handles.DataTypeMenu.Value
               case 1 % 2D-IR Lab 1 & Lab 4 (new)
                   handles = load2DIRlab1(handles);
               case 2 % 2D-IR Lab 4 (old)
                   % Not implemented
               case 3 % Transient 2D-IR
                   handles      = load2DIRlab1(handles);
           end
           guidata(hObject,handles)
           %%% Process the data
           if handles.DataTypeMenu.Value ~= 4
               handles = process2DIR(handles,0);
           end
           % Set t2 delay control values
           handles.Population_delay.String = num2str(handles.t2delays);
           % Ensure that datasets with less t2 delays are loaded correctly
           NewNt2delays = size(handles.t2delays,1);
           if NewNt2delays < OldNt2delays
               handles.Population_delay.Value = 1;
           end
           handles.t2delays = handles.t2delays(:,1);
           % Update handles
           guidata(hObject,handles)
           % Delete progress bar
           delete(handles.WaitBar);
        catch err
           set(handles.ErrorText,'String', string(err.message));
           error=1;
        end
        if error == 0
            % Update the edit controls
            handles = UpdateEdits_2DIR(handles,1);
            % Plot the data in the main window (for preview)
            plot_2DIR(handles,handles.MainAxes);   
            set(handles.MainAxes,'Visible','On')
            % Generate the secondary plot (show phasing data by default)
            SecondaryPlot(handles,handles.SecondaryAxes,'ph');
            set(handles.SecondaryAxes,'Visible','On')
                handles.ShowTimeDomain.Value = 0;
                handles.ShowProbeCalibration.Value = 0;
                handles.ShowPhasing.Value = 1;
            % Enable disabled controls
            EnableControls_2DIR(handles)
            % Update handles
            guidata(hObject,handles)
        end
    end
case 'open'
    % Go to the selected folder
    index_selected = get(handles.DatafoldersList,'Value');
    file_list = get(handles.DatafoldersList,'String');
    datafilename = cellstr(file_list{index_selected}); % Item selected in list box
    % Write datafilename to handles(datafilename) and update
    filenamech = char(datafilename);
    % Go to the selected dir
    cd (filenamech);
    % Update the string
    handles.CurrDir.String = pwd;
    handles.rootdir = handles.CurrDir.String;
    % Put list of subfolders into subFolders
        % Get a list of all files and folders in this folder.
        files = dir(handles.rootdir);
        % Extract only those that are directories.
        subFolders = files([files.isdir]);
        % Don't show the temp folders in the list
        subFolders = subFolders(~contains({subFolders.name},"temp",'IgnoreCase',1));
    % Go to the selected folder
        cd (handles.rootdir);
    % Sort the names and create variables for the list
        dir_struct = subFolders;
        [sorted_names,sorted_index] = sortrows({dir_struct.name}');
        handles.file_names = sorted_names;
        handles.is_dir = [dir_struct.isdir];
        handles.sorted_index = sorted_index;
        guidata(hObject,handles)
    % Update DatafoldersList, removing ".." and "." from the visible list
        handles.file_names = handles.file_names(3:end);
        set(handles.DatafoldersList,'String',handles.file_names,'Value',1)
    guidata(hObject,handles)
end
end


% --- Executes during object creation, after setting all properties.
function DatafoldersList_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ExitButton.
function ExitButton_Callback(hObject, eventdata, handles)
close all


% --- Executes on button press in Refresh.
function Refresh_Callback(hObject, eventdata, handles)
handles.rootdir = handles.CurrDir.String;
% Put list of subfolders into subFolders
    % Get a list of all files and folders in this folder.
    files = dir(handles.rootdir);
    % Extract only those that are directories.
    subFolders = files([files.isdir]);
    % Don't show the temp folders in the list
    subFolders = subFolders(~contains({subFolders.name},"temp",'IgnoreCase',1));
% Go to the folder
    cd (handles.rootdir);
% Sort the names and create variables for the list
    dir_struct = subFolders;
    [sorted_names,sorted_index] = sortrows({dir_struct.name}');
    handles.file_names = sorted_names;
    handles.is_dir = [dir_struct.isdir];
    handles.sorted_index = sorted_index;
    guidata(hObject,handles)
% Update DatafoldersList, removing ".." and "." from the visible list
    handles.file_names = handles.file_names(3:end);
    set(handles.DatafoldersList,'String',handles.file_names,'Value',1)
    guidata(hObject,handles)


% --- Executes on button press in UpDirectory.
function UpDirectory_Callback(hObject, eventdata, handles)
% Go to the folder
cd (handles.CurrDir.String);
% Go one folder Up
cd ..;
% Update the string
handles.CurrDir.String = pwd;
handles.rootdir = handles.CurrDir.String;
% Put list of subfolders into subFolders
    % Get a list of all files and folders in this folder.
    files = dir(handles.rootdir);
    % Extract only those that are directories.
    subFolders = files([files.isdir]);
    % Don't show the temp folders in the list
    subFolders = subFolders(~contains({subFolders.name},"temp",'IgnoreCase',1));
% Go to the selected folder
    cd (handles.rootdir);
% Sort the names and create variables for the list
    dir_struct = subFolders;
    [sorted_names,sorted_index] = sortrows({dir_struct.name}');
    handles.file_names = sorted_names;
    handles.is_dir = [dir_struct.isdir];
    handles.sorted_index = sorted_index;
    guidata(hObject,handles)
% Update DatafoldersList, removing ".." and "." from the visible list
    handles.file_names = handles.file_names(3:end);
    set(handles.DatafoldersList,'String',handles.file_names,'Value',1)
guidata(hObject,handles)


% --- Executes on selection change in DataTypeMenu.
function DataTypeMenu_Callback(hObject, eventdata, handles)
if handles.DataTypeMenu.Value == 3
    handles.Transient2D_text.Visible = 'On';
    handles.Transient2D_mode.Visible = 'On';
else
    handles.Transient2D_text.Visible = 'Off';
    handles.Transient2D_mode.Visible = 'Off';
end

% --- Executes during object creation, after setting all properties.
function DataTypeMenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in ShowTimeDomain.
function ShowTimeDomain_Callback(hObject, eventdata, handles)
td = handles.ShowTimeDomain.Value;
ph = handles.ShowPhasing.Value;
pc = handles.ShowProbeCalibration.Value;

if td==1
    handles.ShowPhasing.Value = 0;
    handles.ShowProbeCalibration.Value = 0;
    SecondaryPlot(handles,handles.SecondaryAxes,'td');
end

if td+ph+pc == 0
    cla(handles.SecondaryAxes,'reset');
    handles.SecondaryAxes.Visible = 'Off';
else
    handles.SecondaryAxes.Visible = 'On';
end

% --- Executes on button press in ShowPhasing.
function ShowPhasing_Callback(hObject, eventdata, handles)
td = handles.ShowTimeDomain.Value;
ph = handles.ShowPhasing.Value;
pc = handles.ShowProbeCalibration.Value;

if ph==1
    handles.ShowTimeDomain.Value = 0;
    handles.ShowProbeCalibration.Value = 0;
    SecondaryPlot(handles,handles.SecondaryAxes,'ph');
end

if td+ph+pc == 0
    cla(handles.SecondaryAxes,'reset');
    handles.SecondaryAxes.Visible = 'Off';
else
    handles.SecondaryAxes.Visible = 'On';
end

% --- Executes on button press in ShowProbeCalibration.
function ShowProbeCalibration_Callback(hObject, eventdata, handles)
td = handles.ShowTimeDomain.Value;
ph = handles.ShowPhasing.Value;
pc = handles.ShowProbeCalibration.Value;

if pc==1
    handles.ShowTimeDomain.Value = 0;
    handles.ShowPhasing.Value = 0;
    SecondaryPlot(handles,handles.SecondaryAxes,'pc');
end

if td+ph+pc == 0
    cla(handles.SecondaryAxes,'reset');
    handles.SecondaryAxes.Visible = 'Off';
else
    handles.SecondaryAxes.Visible = 'On';
end

function PixelNumber_Callback(hObject, eventdata, handles)
if handles.ShowTimeDomain.Value == 1
    SecondaryPlot(handles,handles.SecondaryAxes,'td');
end


% --- Executes during object creation, after setting all properties.
function PixelNumber_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in apodise_method.
function apodise_method_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function apodise_method_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in zeropad_tick.
function zeropad_tick_Callback(hObject, eventdata, handles)


% --- Executes on button press in zeropad_next2k.
function zeropad_next2k_Callback(hObject, eventdata, handles)


function zeropad_factor_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function zeropad_factor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function phase_Npoints_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function phase_Npoints_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in phase_method.
function phase_method_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function phase_method_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function phase_factor_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function phase_factor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BkgSubTick.
function BkgSubTick_Callback(hObject, eventdata, handles)


% --- Executes on button press in ApplyChanges.
function ApplyChanges_Callback(hObject, eventdata, handles)
if handles.DataTypeMenu.Value ~= 4
    % Reprocess the data without reloading
    handles = process2DIR(handles,1);
    guidata(hObject,handles);
    % Clear error string
    set(handles.ErrorText,'String', "");
    % Update the edit controls
    handles = UpdateEdits_2DIR(handles,1);
    % Plot the data in the main window (for preview)
    plot_2DIR(handles,handles.MainAxes);   
    set(handles.MainAxes,'Visible','On')
    % Generate the secondary plot (show phasing data by default)
    SecondaryPlot(handles,handles.SecondaryAxes,'ph');
    set(handles.SecondaryAxes,'Visible','On')
        handles.ShowTimeDomain.Value = 0;
        handles.ShowProbeCalibration.Value = 0;
        handles.ShowPhasing.Value = 1;
else
    handles.WaitBar         = waitbar(0,'Loading data...');
        for j=0:handles.Ndelays
            handles             = load2DIRlab3(handles,j);
            guidata(hObject,handles)
        end
        % Update the edit controls
        handles = UpdateEdits_2DIR(handles,1);
        % Plot the data in the main window (for preview)
        plot_2DIR(handles,handles.MainAxes);   
        set(handles.MainAxes,'Visible','On')

        set(handles.SecondaryAxes,'Visible','Off')
        % Enable disabled controls
        EnableControls_2DIR(handles)
        % Update handles
        guidata(hObject,handles)
end
% Delete progress bar
delete(handles.WaitBar);


% --- Executes on button press in make2Dplot.
function make2Dplot_Callback(hObject, eventdata, handles)
% Copy current plot limits    
    xl= xlim(handles.MainAxes);
    yl= ylim(handles.MainAxes);
% Create figure    
    fh = figure;
    ax = axes(fh);
    ax.Visible='Off';
plot_2DIR(handles,ax);
    fh.Color = [1 1 1];
% Make uniform, consistent format
    ax.FontSize = 12;
    ax.LineWidth = 1;
    ax.TickLength = [0.015 0.035];
    ax.Visible='On';
% Restore axes limits
    xlim(ax,xl)
    ylim(ax,yl)

    
% --- Executes on button press in make3Dplot.
function make3Dplot_Callback(hObject, eventdata, handles)
fh = figure;
ax = axes(fh);
ax.Visible='Off';
SurfPlot_2DIR(handles,ax);
fh.Color = [1 1 1];
% Make uniform, consistent format
ax.FontSize = 12;
ax.LineWidth = 1;
ax.TickLength = [0.015 0.035];
ax.Visible='On';


% --- Executes on button press in plot_Kinetics.
function plot_Kinetics_Callback(hObject, eventdata, handles)
% Get some information from the GUI
plot_pumpdirection  = char(handles.plot_pumpdirection.String{handles.plot_pumpdirection.Value});
ProbeAxis           = handles.ProbeAxis;
PumpAxis            = handles.PumpAxis;
Ndelays             = handles.Ndelays;
t2delays            = handles.t2delays;
PROC_2D_DATA        = handles.PROC_2D_DATA;

% Get the desired values to plot from user
    handles.SelTraces = [];
    switch handles.InteractiveModeTick.Value
        case 1
            % Select points interactively, to finish hit RETURN
            handles = SelectTraces(handles,0);
            if isempty(handles.SelTraces)
                return
            end
            EnT = 0;
        case 0
            handles.SelTraces = inputdlg('Enter the coordinates of the points to plot along t2 (format: x1,y1;x2,y2...):',...
                 'Input desired coordinates', [1 60]);
            if isempty(handles.SelTraces)
                return
            end
            % Plot special cases
            if strcmp(handles.SelTraces,'EnT-GSB')
                handles.SelTraces = [1979,1979;1979,2028;2028,2028;2028,1979];
                EnT = 1;
            elseif strcmp(handles.SelTraces,'EnT-ESA')
                handles.SelTraces = [1979,1967;1979,2017;2028,2017;2028,1967];
                EnT = 1;
            elseif strcmp(handles.SelTraces,'Andrea')
                handles.SelTraces = [1985,1985;1985,2062;2062,2062;2062,1985];
                EnT = 2;
            else
                handles.SelTraces = str2num(handles.SelTraces{:});
                EnT = 0;
            end
    end
    L = size(handles.SelTraces,1);
    
% Separate the values into pump and probe, according to how the data is plotted    
switch plot_pumpdirection
    case 'Horizontal'
        pump_search     = handles.SelTraces(:,1); 
        probe_search    = handles.SelTraces(:,2);
    case 'Vertical'
        pump_search     = handles.SelTraces(:,2); 
        probe_search    = handles.SelTraces(:,1);
end

% Find the desired values in the pump and probe axis
for i=1:L
    probe_index(i)          = findClosestId2Val(ProbeAxis,probe_search(i));
    for m=1:Ndelays
        pump_index(m,i)     = findClosestId2Val(PumpAxis{m,1},pump_search(i));
    end
end

% Initialise variables
caption = cell(L,1);
kindata = zeros(Ndelays,L);
    
for i=1:L
    % Get the values from the 2D arrays for each population delay and plot them
    % The data is by default in the format (w1,w3) (according to process2DIR.m)
    % The vector ydata will have the following format: (assuming N traces and M t2 delays)
    %         t2 delay 1: [Z1 Z2 Z3 ... Zn]
    %         t2 delay 2: [Z1 Z2 Z3 ... Zn]
    %         t2 delay 3: [Z1 Z2 Z3 ... Zn]
    %             ...
    %         t2 delay m: [Z1 Z2 Z3 ... Zn]
    
    % Get the data and normalise (if needed)
    for m=1:Ndelays
        kindata(m,i) = PROC_2D_DATA{m,1}(pump_index(m,i),probe_index(i));
    end
    
    switch handles.Normalise.Value
        case 0
            label = '2D signal (a.u.)';
            caption{i} = ['(' num2str(round(PumpAxis{1,1}(pump_index(1,i)))) ', ' num2str(ProbeAxis(probe_index(i))) ') cm^{-1}'];
        case 1
            t0_index = findClosestId2Val(t2delays,0);
            maxval = max(kindata(t0_index:end,i));
            minval = min(kindata(t0_index:end,i));
            if abs(maxval) >= abs(minval)
                kindata(:,i) = kindata(:,i)./maxval;
                caption{i} = ['(' num2str(round(PumpAxis{1,1}(pump_index(1,i)))) ', ' num2str(ProbeAxis(probe_index(i))) ') cm^{-1}'];
            else
                kindata(:,i) = kindata(:,i)./minval;
                caption{i} = ['(' num2str(round(PumpAxis{1,1}(pump_index(1,i)))) ', ' num2str(ProbeAxis(probe_index(i))) ') cm^{-1} \times -1'];
            end
            label = 'Normalised 2D signal (a.u.)';
    end
end



% Prepare plot for EnT
if EnT ~= 0
    if abs(max(kindata(:))) <= abs(min(kindata(:)))
        kindata = -kindata;
    end
    
    diagFW          = kindata(:,1)./max(max(abs(kindata(:,1:2))));
    diagBW          = kindata(:,3)./max(max(abs(kindata(:,3:4))));
    xpeakFW         = 10.*kindata(:,2)./max(max(abs(kindata(:,1:2))));
    xpeakBW         = 10.*kindata(:,4)./max(max(abs(kindata(:,3:4))));
    time            = handles.t2delays;
    % Create figure
    fh              = figure;
    fh.Units        = 'normalized';
    fh.Position(2)  = 0.1;
    fh.Position(4)  = fh.Position(4)*2;

    fh.Color        = [1 1 1];
    fh.Units        = 'pixels';
    ax_FW           = subplot(2,1,1);
    ax_BW           = subplot(2,1,2);
    ax_FW.FontSize  = 16;
    ax_BW.FontSize  = 16;

    box(ax_FW,'on');
    box(ax_BW,'on');

    hold(ax_FW,'on');
    hold(ax_BW,'on');
    switch EnT
        case 1
            diag1 = 'Diagonal Re(^{13}CO)';
            diag2 = 'Diagonal Re(^{12}CO)';
            xpeak1= 'Re(^{13}CO) \rightarrow Re(^{12}CO) \rm{\times10}';
            xpeak2= 'Re(^{13}CO) \leftarrow Re(^{12}CO) \rm{\times10}';
            titFW = 'Forward energy transfer';
            titBW = 'Backward energy transfer';
        case 2
            diag2 = 'Diagonal SCN^{-}';
            diag1 = 'Diagonal S^{13}C^{15}N^{-}';
            xpeak2= 'SCN^{-} \rightarrow S^{13}C^{15}N^{-} \rm{\times 10}';
            xpeak1= 'SCN^{-} \leftarrow S^{13}C^{15}N^{-} \rm{\times 10}';
            titBW = 'Downhill transfer';
            titFW = 'Uphill transfer';
    end
    % Plot the diagonal peaks
    plot(ax_FW,time,diagFW,'-or','LineWidth',1,'DisplayName',diag1);
    plot(ax_BW,time,diagBW,'-ob','LineWidth',1,'DisplayName',diag2);

    % Plot the cross peaks
    plot(ax_FW,time,xpeakFW,'-^r','LineWidth',1,'DisplayName',xpeak1);
    plot(ax_BW,time,xpeakBW,'-vb','LineWidth',1,'DisplayName',xpeak2);

    % Set axes limits
    axis(ax_FW,'tight');
    axis(ax_BW,'tight');
    
    %%% Nice formatting
    % Titles
    title(ax_FW,titFW,'FontSize',16);
    title(ax_BW,titBW,'FontSize',16);

    % Axis labels
    xlabel(ax_FW,'t_2 delay (ps)','FontWeight','bold','FontSize',16);
    xlabel(ax_BW,'t_2 delay (ps)','FontWeight','bold','FontSize',16);

    ylabel(ax_FW,'Normalised 2D signal','FontWeight','bold','FontSize',16);
    ylabel(ax_BW,'Normalised 2D signal','FontWeight','bold','FontSize',16);

    % Axis limits
    xlim(ax_FW,[0 time(end)]);
    xlim(ax_BW,[0 time(end)]);

    ylim(ax_FW,[-0.1 1.1]);
    ylim(ax_BW,[-0.1 1.1]);
    % Add zero line
    hline_FW = refline(ax_FW,[0 0]); hline_FW.Color = [0.5 0.5 0.5];
    set(get(get(hline_FW,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    hline_BW = refline(ax_BW,[0 0]); hline_BW.Color = [0.5 0.5 0.5];
    set(get(get(hline_BW,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    
    % Legends
    lh_FW = legend(ax_FW,'show');
    legend(ax_FW,'boxoff','FontWeight','bold')
    legend(ax_FW,{},'FontWeight','bold')

    lh_BW = legend(ax_BW,'show');
    legend(ax_BW,'boxoff')
    legend(ax_BW,{},'FontWeight','bold')
    % Make figure resizable
    fh.Units        = 'normalized';
else
    % Create a new figure with consistent format
    fh = figure();
    fh.Position(3)  = 800;
    fh.Position(4)  = 425;
    fh.Color        = [1 1 1];
    % Define the axes
    axes2 = axes('Parent',fh);
    axes(axes2);
    cmap=colormap(othercolor('Mrainbow',L));
    % Plot the data
    for n=1:L
       plot(axes2,handles.t2delays,kindata(:,n),'LineWidth',2,'Marker','o','MarkerSize',2,'color',cmap(n,:));
       hold on
    end
    %%% Nice formatting
    set(gca,'FontSize',14)
    xlabel('t_{2} delay (ps)','FontSize',14,'FontWeight','bold');
    ylabel(label,'FontSize',14,'FontWeight','bold')
    % title([handles.datafilename;'WAITING TIME KINETICS';''],'Interpreter','none','FontSize',10)
    axis tight

    % Show only positive t2 times
    xlim([0 max(t2delays)]);

    % Create legend
    legend(gca,caption,'FontSize',12)
    legend('boxoff')
    legend('Location','northeast')

    % Set linear or log X scale
    set(axes2,'xscale','lin');

    % Add zero line
    hline = refline(axes2,[0 0]); hline.Color = [0.5 0.5 0.5];
    set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

    axes2.Units     = 'pixels';
    axes2.Position  = [75 75 675 320];
    axes2.Units     = 'normalized';
end


% Decide whether to save the plotted traces or not
% if handles.DoSaveTraces==1
%     wavenumbers=transpose(handles.cmprobe(k));
    filename    = char(strcat(handles.CurrDir.String,filesep,handles.datafilename,'_traces.dat'));
%     data        = [[[0;PumpAxis{1,1}(pump_index(1,:))],[0;ProbeAxis(probe_index)]]';[handles.t2delays(2:end),kindata(2:end,:)]];
    data        = [[[0;PumpAxis{1,1}(pump_index(1,:))],[0;ProbeAxis(probe_index)]]';[handles.t2delays,[diagFW xpeakFW diagBW xpeakBW]]];
    dlmwrite(filename,data);
% end
guidata(hObject,handles)


% --- Executes on button press in plot_Slices.
function plot_Slices_Callback(hObject, eventdata, handles)
% Ask the user what to plot
slice_options = {'Diagonal','Along fixed pump WL','Along fixed probe WL','Integrate along pump axis','Integrate along probe axis','Along several pump WL','Along several probe WL'};
[slice_typeindx,doplot] = listdlg('ListString',slice_options,'OKstring','Plot','SelectionMode','single','ListSize',[150,120],'PromptString','Select slice type to plot:');

if doplot == 0
    return
end

% Get some information from the GUI
plot_pumpdirection  = char(handles.plot_pumpdirection.String{handles.plot_pumpdirection.Value});
ProbeAxis           = handles.ProbeAxis;
PumpAxis            = handles.PumpAxis;
Ndelays             = handles.Ndelays;
t2delays            = handles.t2delays;
PROC_2D_DATA        = handles.PROC_2D_DATA;
PopDelay            = handles.Population_delay.Value;
Normalise           = handles.Normalise.Value;

switch slice_options{slice_typeindx}
    case 'Diagonal'
        % Initialise variables
        pump_indexes    = cell(Ndelays,1);
        data            = zeros(length(ProbeAxis),Ndelays);
        
        % Get the data to be plotted
        for m=1:Ndelays
            % Get the indices of the probe wavelengths in the pump axis
            pump_indexes{m}     = findClosestId2Val(PumpAxis{m,1},transpose(ProbeAxis));
            probe_indexes       = 1:length(ProbeAxis);
            % Get the data
            data(:,m)           = diag(PROC_2D_DATA{m,1}(pump_indexes{m},probe_indexes));
        end
        
        % Create a new figure
        fh = figure();

        % Define the axes
        handles.axes2 = axes(fh);
        
        % Plot the data
        cmap=colormap(othercolor('Mrainbow',Ndelays));
        for m=1:Ndelays
           plot(handles.axes2,ProbeAxis,data(:,m),'-','LineWidth',2,'MarkerSize',2,'color',cmap(m,:),'DisplayName',['t_{2} = ' num2str(t2delays(m),'%.3g') ' ps']);
           hold on
        end
        
        % Save plot details
        plot_title  = 'DIAGONAL CUTS';
        x_axis      = 'Probe';
        
        %%% Nice formatting
        % Size and colour
        fh.Color = [1 1 1];
        currsize = fh.Position;
        fh.Position = [currsize(1:2) 900 425];
        handles.axes2.Units = 'Pixels';
        currsize = handles.axes2.OuterPosition;
        handles.axes2.OuterPosition = [currsize(1:2) 950 420];
        handles.axes2.Units = 'Normalized';

        % Title, axis labels and legend
        handles.axes2.FontSize = 14;
        xlabel(handles.axes2,[x_axis ' wavelength (cm^{-1})'],'FontSize',14,'FontWeight','bold');
        ylabel(handles.axes2,'2D signal (a.u.)','FontSize',14,'FontWeight','bold')
        title(handles.axes2,[handles.datafilename;plot_title;''],'Interpreter','none','FontSize',10)
        axis tight

        % Create legend
        legend(gca);
        legend('boxoff')
        legend('Location','eastoutside')

        % Add zero line
        hline = refline(handles.axes2,[0 0]);
        hline.Color = [0.5 0.5 0.5];
        set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        
    case 'Along fixed pump WL'
        switch handles.InteractiveModeTick.Value
            case 1
                % Select points interactively, to finish hit RETURN
                handles = SelectTraces(handles,0); 
                % Separate the values into pump and probe, according to how the data is plotted    
                switch plot_pumpdirection
                    case 'Horizontal'
                        pump_search    = handles.SelTraces(:,1); 
                    case 'Vertical'
                        pump_search    = handles.SelTraces(:,2); 
                end
            case 0
                handles.SelTraces = inputdlg('Enter the pump wavenumbers to plot:',...
                     'Input pump wavenumbers to plot:', [1 60]);
                pump_search = str2num(handles.SelTraces{:})';
        end
        L = length(pump_search);

        % Initialise variables
        pump_indexes    = cell(Ndelays,1);
        data            = zeros(length(ProbeAxis),Ndelays,L);
        
        % Get the data to be plotted
        for m=1:Ndelays
            % Get the indices of the probe wavelengths
            pump_indexes{m}     = findClosestId2Val(PumpAxis{m,1},pump_search);
            probe_indexes       = 1:length(ProbeAxis);
            % Get the data
            for p=1:L
                data(:,m,p)     = transpose(PROC_2D_DATA{m,1}(pump_indexes{m}(p),probe_indexes));
            end
        end
        
        % Plot each pump wavelength in one new figure
        for p=1:L
            % Create a new figure
            fh = figure();

            % Define the axes
            handles.axes2 = axes(fh);

            % Plot the data
            cmap=colormap(othercolor('Mrainbow',Ndelays));
            for m=1:Ndelays
               plot(handles.axes2,ProbeAxis,data(:,m,p),'-','LineWidth',2,'MarkerSize',2,'color',cmap(m,:),'DisplayName',['t_{2} = ' num2str(t2delays(m),'%.3g') ' ps']);
               hold on
            end

            % Save plot details
            plot_title  = ['TRANSIENTS AT ' num2str(PumpAxis{1,1}(pump_indexes{1}(p)),'%.4g') ' cm^{-1} (pump)'];
            x_axis      = 'Probe';
            
            %%% Nice formatting
            % Size and colour
            fh.Color = [1 1 1];
            currsize = fh.Position;
            fh.Position = [currsize(1:2) 900 425];
            handles.axes2.Units = 'Pixels';
            currsize = handles.axes2.OuterPosition;
            handles.axes2.OuterPosition = [currsize(1:2) 950 420];
            handles.axes2.Units = 'Normalized';

            % Title, axis labels and legend
            handles.axes2.FontSize = 14;
            xlabel(handles.axes2,[x_axis ' wavelength (cm^{-1})'],'FontSize',14,'FontWeight','bold');
            ylabel(handles.axes2,'2D signal (a.u.)','FontSize',14,'FontWeight','bold')
            title(handles.axes2,plot_title,'FontSize',10)
            axis tight

            % Create legend
            legend(gca);
            legend('boxoff')
            legend('Location','eastoutside')

            % Add zero line
            hline = refline(handles.axes2,[0 0]);
            hline.Color = [0.5 0.5 0.5];
            set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        end
    
    case 'Along fixed probe WL'
        % Get the desired values to plot from user
        handles.SelTraces = [];
        switch handles.InteractiveModeTick.Value
            case 1
                % Select points interactively, to finish hit RETURN
                handles = SelectTraces(handles,0); 
                % Separate the values into pump and probe, according to how the data is plotted    
                switch plot_pumpdirection
                    case 'Horizontal'
                        probe_search    = handles.SelTraces(:,2); 
                    case 'Vertical'
                        probe_search    = handles.SelTraces(:,1); 
                end
            case 0
                handles.SelTraces = inputdlg('Enter the probe wavenumbers to plot:',...
                     'Input probe wavenumbers to plot:', [1 60]);
                probe_search = str2num(handles.SelTraces{:});
        end
        L = length(probe_search);
        
        % Get the segment of the pump axis
        for m=1:Ndelays
            if handles.CutPlot_tick.Value
                min_pump        = findClosestId2Val(PumpAxis{m,1},min(ProbeAxis));
                max_pump        = findClosestId2Val(PumpAxis{m,1},max(ProbeAxis));
                pump_indexes{m} = min_pump:1:max_pump;
            else
                pump_indexes{m} = 1:1:length(PumpAxis{m,1});
            end
        end
        
        % Initialise variables
        data                    = zeros(length(PumpAxis{1,1}(pump_indexes{m})),Ndelays,L);
              
        % Get the probe indexes
        probe_indexes           = findClosestId2Val(ProbeAxis,probe_search);
        
        % Get the data to be plotted
        for m=1:Ndelays
            for p=1:L
                data(:,m,p)     = PROC_2D_DATA{m,1}(pump_indexes{m},probe_indexes(p));
            end
        end
        
        % Plot each probe wavelength in one new figure
        for p=1:L
            % Create a new figure
            fh = figure();

            % Define the axes
            handles.axes2 = axes(fh);

            % Plot the data
            cmap=colormap(othercolor('Mrainbow',Ndelays));
            for m=1:Ndelays
               plot(handles.axes2,PumpAxis{m,1}(pump_indexes{m}),data(:,m,p),'-','LineWidth',2,'MarkerSize',2,'color',cmap(m,:),'DisplayName',['t_{2} = ' num2str(t2delays(m),'%.3g') ' ps']);
               hold on
            end

            % Save plot details
            plot_title  = ['TRANSIENTS AT ' num2str(ProbeAxis(probe_indexes(p)),'%.4g') ' cm^{-1} (probe)'];
            x_axis      = 'Pump';
            
            %%% Nice formatting
            % Size and colour
            fh.Color = [1 1 1];
            currsize = fh.Position;
            fh.Position = [currsize(1:2) 900 425];
            handles.axes2.Units = 'Pixels';
            currsize = handles.axes2.OuterPosition;
            handles.axes2.OuterPosition = [currsize(1:2) 950 420];
            handles.axes2.Units = 'Normalized';

            % Title, axis labels and legend
            handles.axes2.FontSize = 14;
            xlabel(handles.axes2,[x_axis ' wavelength (cm^{-1})'],'FontSize',12,'FontWeight','bold');
            ylabel(handles.axes2,'2D signal (a.u.)','FontSize',12,'FontWeight','bold')
            title(handles.axes2,plot_title,'FontSize',10)
            axis tight

            % Create legend
            legend(gca);
            legend('boxoff')
            legend('Location','eastoutside')

            % Add zero line
            hline = refline(handles.axes2,[0 0]);
            hline.Color = [0.5 0.5 0.5];
            set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        end

    case 'Integrate along pump axis'
        switch handles.InteractiveModeTick.Value
            case 0
                % Ask user for pump range (default is same as probe range)
                defaults  = {num2str([ProbeAxis(1) ProbeAxis(end)])};
                SelTraces = inputdlg('Enter the pump range to integrate:',...
                             'Input pump wavenumber range to integrate:', [1 60],defaults);
                SelTraces = str2num(SelTraces{:});
            case 1
                SelTraces = [ProbeAxis(1) ProbeAxis(end)];
        end
        
        % Get the indices of the pump wavelengths to integrate
        pump_indexes    = findClosestId2Val(PumpAxis{1,1},SelTraces);
        % Initialize variables
        data            = zeros(length(ProbeAxis),Ndelays);
        % Get the data to be plotted
        for m=1:Ndelays
            data(:,m)   = transpose(sum(PROC_2D_DATA{m,1}(min(pump_indexes):max(pump_indexes),:),1));
        end
        
        % Create a new figure
        fh = figure();

        % Define the axes
        handles.axes2 = axes(fh);

        % Plot the data
        cmap=colormap(othercolor('Mrainbow',Ndelays));
        for m=1:Ndelays
           plot(handles.axes2,ProbeAxis,data(:,m),'-','LineWidth',2,'MarkerSize',2,'color',cmap(m,:),'DisplayName',[num2str(t2delays(m),'%.3g') ' ps']);
           hold on
        end

        % Save plot details
        plot_title  = ['INTEGRATED 2D SIGNAL FOR PUMP WL' num2str(PumpAxis{1,1}(min(pump_indexes))) ' to ' num2str(PumpAxis{1,1}(max(pump_indexes))) ' cm^{-1}'];
        x_axis      = 'Probe';

        %%% Nice formatting
        % Size and colour
        fh.Color = [1 1 1];
        currsize = fh.Position;
        fh.Position = [currsize(1:2) 900 425];
        handles.axes2.Units = 'Pixels';
        currsize = handles.axes2.OuterPosition;
        handles.axes2.OuterPosition = [currsize(1:2) 950 420];
        handles.axes2.Units = 'Normalized';

        % Title, axis labels and legend
        handles.axes2.FontSize = 14;
        xlabel(handles.axes2,[x_axis ' wavelength (cm^{-1})'],'FontSize',12,'FontWeight','bold');
        ylabel(handles.axes2,'2D signal (a.u.)','FontSize',12,'FontWeight','bold')
        title(handles.axes2,plot_title,'FontSize',10)
        axis tight

        % Create legend
        legend(gca);
        legend('boxoff')
        legend('Location','eastoutside')

        % Add zero line
        hline = refline(handles.axes2,[0 0]);
        hline.Color = [0.5 0.5 0.5];
        set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
            
    case 'Integrate along probe axis'
        switch handles.InteractiveModeTick.Value
            case 0
                % Ask user for pump range (default is same as probe range)
                defaults          = {num2str([ProbeAxis(1) ProbeAxis(end)])};
                SelTraces = inputdlg('Enter the pump range to show:',...
                             'Input pump wavenumber range to show:', [1 60],defaults);
                SelTraces = str2num(SelTraces{:});
            case 1
                SelTraces = [ProbeAxis(1) ProbeAxis(end)];
        end
        
        % Get the indices of the pump wavelengths
        pump_indexes    = findClosestId2Val(PumpAxis{1,1},SelTraces);
        % Initialize variables
        data            = zeros((max(pump_indexes)-min(pump_indexes)+1),Ndelays);
        % Get the data to be plotted
        for m=1:Ndelays
            data(:,m)   = sum(PROC_2D_DATA{m,1}(min(pump_indexes):max(pump_indexes),:),2);
        end
        
        % Create a new figure
        fh = figure();

        % Define the axes
        handles.axes2 = axes(fh);

        % Plot the data
        cmap=colormap(othercolor('Mrainbow',Ndelays));
        for m=1:Ndelays
           plot(handles.axes2,PumpAxis{1,1}(min(pump_indexes):max(pump_indexes)),data(:,m),'-','LineWidth',2,'MarkerSize',2,'color',cmap(m,:),'DisplayName',['t_{2} = ' num2str(t2delays(m),'%.3g') ' ps']);
           hold on
        end

        % Save plot details
        plot_title  = ['2D SIGNAL, INTEGRAL ACROSS ALL PROBE WAVELENGTHS'];
        x_axis      = 'Pump';

        %%% Nice formatting
        % Size and colour
        fh.Color = [1 1 1];
        currsize = fh.Position;
        fh.Position = [currsize(1:2) 900 425];
        handles.axes2.Units = 'Pixels';
        currsize = handles.axes2.OuterPosition;
        handles.axes2.OuterPosition = [currsize(1:2) 950 420];
        handles.axes2.Units = 'Normalized';

        % Title, axis labels and legend
        handles.axes2.FontSize = 14;
        xlabel(handles.axes2,[x_axis ' wavelength (cm^{-1})'],'FontSize',12,'FontWeight','bold');
        ylabel(handles.axes2,'2D signal (a.u.)','FontSize',12,'FontWeight','bold')
        title(handles.axes2,plot_title,'FontSize',10)
        axis tight

        % Create legend
        legend(gca);
        legend('boxoff')
        legend('Location','eastoutside')

        % Add zero line
        hline = refline(handles.axes2,[0 0]);
        hline.Color = [0.5 0.5 0.5];
        set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    case 'Along several pump WL'
        switch handles.InteractiveModeTick.Value
            case 1
                % Select points interactively, to finish hit RETURN
                handles = SelectTraces(handles,0); 
                % Separate the values into pump and probe, according to how the data is plotted    
                switch plot_pumpdirection
                    case 'Horizontal'
                        pump_search    = handles.SelTraces(:,1); 
                    case 'Vertical'
                        pump_search    = handles.SelTraces(:,2); 
                end
            case 0
                handles.SelTraces = inputdlg('Enter the pump wavenumbers to plot:',...
                     'Input pump wavenumbers to plot:', [1 60]);
                pump_search = str2num(handles.SelTraces{:})';
        end
        
        % Get the data to be plotted
        % Get the indices of the probe wavelengths
        pump_indexes        = unique(findClosestId2Val(PumpAxis{PopDelay,1},pump_search));
        probe_indexes       = 1:length(ProbeAxis);
        L                   = length(pump_indexes);
        % Initialise variables
        data            = zeros(length(ProbeAxis),L);
        % Get the data
        for p=1:L
            data(:,p)       = transpose(PROC_2D_DATA{PopDelay,1}(pump_indexes(p),probe_indexes));
        end
        
        if Normalise
            data            = data./max(abs(data),[],1);
        end
        
        % Create a new figure
        fh = figure();

        % Define the axes
        handles.axes2 = axes(fh);

        % Plot the data
        cmap=colormap(othercolor('Mrainbow',L));
        for p=1:L
           plot(handles.axes2,ProbeAxis,data(:,p),'-','LineWidth',2,'MarkerSize',2,'color',cmap(p,:),'DisplayName',[num2str(PumpAxis{PopDelay,1}(pump_indexes(p)),'%.5g') ' cm^{-1}']);
           hold on
        end

        % Save plot details
        plot_title  = ['CUTS ACROSS PUMP WAVELENGTHS AT t_2 = ' num2str(t2delays(PopDelay),'%.3g') ' ps'];
        x_axis      = 'Probe';

        %%% Nice formatting
        % Size and colour
        fh.Color = [1 1 1];
        currsize = fh.Position;
        fh.Position = [currsize(1:2) 900 425];
        handles.axes2.Units = 'Pixels';
        currsize = handles.axes2.OuterPosition;
        handles.axes2.OuterPosition = [currsize(1:2) 950 420];
        handles.axes2.Units = 'Normalized';

        % Title, axis labels and legend
        handles.axes2.FontSize = 14;
        xlabel(handles.axes2,[x_axis ' wavelength (cm^{-1})'],'FontSize',14,'FontWeight','bold');
        ylabel(handles.axes2,'2D signal (a.u.)','FontSize',14,'FontWeight','bold')
        title(handles.axes2,plot_title,'FontSize',10)
        axis tight

        % Create legend
        legend(gca);
        legend('boxoff')
        legend('Location','eastoutside')

        % Add zero line
        hline = refline(handles.axes2,[0 0]);
        hline.Color = [0.5 0.5 0.5];
        set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
end


% --- Executes on button press in PowerDependence.
function PowerDependence_Callback(hObject, eventdata, handles)

% --- Executes on button press in Fit.
function Fit_Callback(hObject, eventdata, handles)
handles = Gaussian2D_analysis(handles);
guidata(hObject, handles);

% --- Executes on button press in ShiftT2.
function ShiftT2_Callback(hObject, eventdata, handles)

% --- Executes on button press in SpectralDiffusion.
function SpectralDiffusion_Callback(hObject, eventdata, handles)
handles = SpectralDiffusion_analysis(handles);
plot_2DIR(handles,handles.MainAxes);
guidata(hObject, handles);

% --- Executes on button press in IntegralDynamics.
function IntegralDynamics_Callback(hObject, eventdata, handles)
% Ask the user what to plot
integral_options    = {'Peak volume kinetics','Peak max-min kinetics','Peak volume difference kinetics','Noise analysis'};
[integral_typeindx,doplot] = listdlg('ListString',integral_options,'OKstring','Plot','SelectionMode','single','ListSize',[200,100],'PromptString','Select kinetic analysis type:');

if doplot == 0
    return
end

% Get some information from the GUI
plot_pumpdirection  = char(handles.plot_pumpdirection.String{handles.plot_pumpdirection.Value});
ProbeAxis           = handles.ProbeAxis;
PumpAxis            = handles.PumpAxis;
Ndelays             = handles.Ndelays;
t2delays            = handles.t2delays;
PROC_2D_DATA        = handles.PROC_2D_DATA;
Axes                = handles.MainAxes;

switch integral_options{integral_typeindx}
    case 'Peak volume kinetics'
        % Get the desired values to plot from user
        handles.SelTraces = [];
        
        switch handles.InteractiveModeTick.Value
            case 1
                % Select points interactively, to finish hit RETURN
                [RegionPositions,L] = SelectRegions(Axes);
                % Separate the values into pump and probe, according to how the data is plotted
                % i.e. convert from [x_min x_max y_min y_max] to [pump_min pump_max probe_min probe_max]
                switch plot_pumpdirection
                    case 'Horizontal'
                        pump_ranges     = RegionPositions(:,1:2); 
                        probe_ranges    = RegionPositions(:,3:4);
                    case 'Vertical'
                        pump_ranges     = RegionPositions(:,3:4); 
                        probe_ranges    = RegionPositions(:,1:2);
                end
            case 0
                RegionPositions_cell    = inputdlg({'Enter the pump regions (Format: min1 max1; min2 max2; ...): ','Enter the probe regions (Format: min1 max1; min2 max2; ...):'},'Define regions to integrate', [1 80]);
                pump_ranges              = str2num(RegionPositions_cell{1});
                probe_ranges             = str2num(RegionPositions_cell{2});
                L = size(pump_ranges,1);
                M = size(probe_ranges,1);
                if L ~= M
                    warndlg('Number of pump and probe ranges is different!')
                    return
                end
        end

        % Initialise variables
        pump_idx    = zeros(L,2);
        probe_idx   = zeros(L,2);
        kindata     = zeros(Ndelays,L);
        
        % Get the data to be plotted
        for m=1:Ndelays
            % Get the values from the 2D arrays for each population delay on the selected regions
            % The data is by default in the format (w1,w3) (according to process2DIR.m)
            for p=1:L
                % Get the indices of the probe wavelengths
                pump_idx(p,:)   = findClosestId2Val(PumpAxis{m,1},pump_ranges(p,:));
                probe_idx(p,:)  = findClosestId2Val(ProbeAxis,probe_ranges(p,:));
                tempdata        = PROC_2D_DATA{m,1}(pump_idx(p,1):pump_idx(p,2),probe_idx(p,1):probe_idx(p,2));
                kindata(m,p)    = sum(tempdata(:));
            end
        end
        % Normalise the data
        if handles.Normalise.Value == 1
            kindata    = kindata./max(abs(kindata),[],1);
            label      = 'Normalised Peak volume (a.u.)';
        else
            label      = 'Peak volume (a.u.)';
        end
        
        % Create a new figure
        fh          = figure();
        fh.Color    = [1 1 1];

        % Define the axes
        axes2       = axes('Parent',fh);

        % Plot the data
        cmap=colormap(othercolor('Mrainbow',L));
        for n=1:L
           plotname = ['(' num2str(PumpAxis{m,1}(pump_idx(n,1)),'%.4g') '-' num2str(PumpAxis{m,1}(pump_idx(n,2)),'%.4g') '; ' num2str(ProbeAxis(probe_idx(n,1)),'%.4g') '-' num2str(ProbeAxis(probe_idx(n,2)),'%.4g') ') cm^{-1}'];
           plot(axes2,t2delays,kindata(:,n),'LineWidth',2,'Marker','o','MarkerSize',2,'color',cmap(n,:),'DisplayName',plotname);
           hold on
        end

        %%% Nice formatting
        set(gca,'FontSize',14)
        xlabel('t_{2} delay (ps)','FontSize',14,'FontWeight','bold');
        ylabel(label,'FontSize',14,'FontWeight','bold')
        title([handles.datafilename;'WAITING TIME KINETICS';''],'Interpreter','none','FontSize',10)
        axis tight

        % Show only positive t2 times
        xlim([0.25 max(t2delays)]);

        % Create legend
        legend(gca,{},'FontSize',12)
        legend('boxoff')
        legend('Location','northeast')

        % Set linear or log X scale
        set(axes2,'xscale','lin');

        % Add zero line
        hline = refline(axes2,[0 0]); hline.Color = [0.5 0.5 0.5];
        set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

        % % Decide whether to save the plotted traces or not
        % if handles.DoSaveTraces==1
        %     wavenumbers=transpose(handles.cmprobe(k));
        %     filename=char(strcat(handles.CurrDir.String,filesep,handles.datafilename,'_traces.dat'));
        %     dlmwrite(filename,[[0;handles.delays] [wavenumbers;kindata]])
        % end
        guidata(hObject,handles)
    case 'Noise analysis'
        % Get the desired values to plot from user
        handles.SelTraces = [];
        
        switch handles.InteractiveModeTick.Value
            case 1
                % Select points interactively, to finish hit RETURN
                [RegionPositions,L] = SelectRegions(Axes);
                % Separate the values into pump and probe, according to how the data is plotted
                % i.e. convert from [x_min x_max y_min y_max] to [pump_min pump_max probe_min probe_max]
                switch plot_pumpdirection
                    case 'Horizontal'
                        pump_ranges     = RegionPositions(:,1:2); 
                        probe_ranges    = RegionPositions(:,3:4);
                    case 'Vertical'
                        pump_ranges     = RegionPositions(:,3:4); 
                        probe_ranges    = RegionPositions(:,1:2);
                end
            case 0
                RegionPositions_cell    = inputdlg({'Enter the pump regions (Format: min1 max1; min2 max2; ...): ','Enter the probe regions (Format: min1 max1; min2 max2; ...):'},'Define regions to integrate', [1 80]);
                pump_ranges              = str2num(RegionPositions_cell{1});
                probe_ranges             = str2num(RegionPositions_cell{2});
                L = size(pump_ranges,1);
                M = size(probe_ranges,1);
                if L ~= M
                    warndlg('Number of pump and probe ranges is different!')
                    return
                end
        end

        % Initialise variables
        pump_idx    = zeros(L,2);
        probe_idx   = zeros(L,2);
        noise       = zeros(Ndelays,L);
        
        % Get the data to be plotted
        for m=1:Ndelays
            % Get the values from the 2D arrays for each population delay on the selected regions
            % The data is by default in the format (w1,w3) (according to process2DIR.m)
            for p=1:L
                % Get the indices of the probe wavelengths
                pump_idx(p,:)   = findClosestId2Val(PumpAxis{m,1},pump_ranges(p,:));
                probe_idx(p,:)  = findClosestId2Val(ProbeAxis,probe_ranges(p,:));
                tempdata        = PROC_2D_DATA{m,1}(pump_idx(p,1):pump_idx(p,2),probe_idx(p,1):probe_idx(p,2));
                noise(m,p)      = std(tempdata(:));
            end
        end
        warndlg(['The noise level is: ' num2str(noise)]);
        guidata(hObject,handles)        
end


% --- Executes on button press in InteractiveModeTick.
function InteractiveModeTick_Callback(hObject, eventdata, handles)


% --- Executes on button press in Normalise.
function Normalise_Callback(hObject, eventdata, handles)


% --- Executes on button press in Probe_Calibration.
function Probe_Calibration_Callback(hObject, eventdata, handles)


% --- Executes on button press in PumpCorrection_tick.
function PumpCorrection_tick_Callback(hObject, eventdata, handles)


function editWLmin_Callback(hObject, eventdata, handles)
maxWL       = str2double(handles.editWLmax.String);
minWL       = str2double(handles.editWLmin.String);
plotaxis    = handles.MainAxes;
plot_pumpdirection  = char(handles.plot_pumpdirection.String{handles.plot_pumpdirection.Value});
if handles.EditProbeAxis_tick.Value == 1
    switch plot_pumpdirection
    case 'Vertical' 
        xlim(plotaxis,[minWL,maxWL]);
    case 'Horizontal'
        ylim(plotaxis,[minWL,maxWL]);
    end
else
    switch plot_pumpdirection
    case 'Vertical' 
        ylim(plotaxis,[minWL,maxWL]);
    case 'Horizontal'
        xlim(plotaxis,[minWL,maxWL]);
    end
end


% --- Executes during object creation, after setting all properties.
function editWLmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editWLmax_Callback(hObject, eventdata, handles)
maxWL       = str2double(handles.editWLmax.String);
minWL       = str2double(handles.editWLmin.String);
plotaxis    = handles.MainAxes;
plot_pumpdirection  = char(handles.plot_pumpdirection.String{handles.plot_pumpdirection.Value});
if handles.EditProbeAxis_tick.Value == 1
    switch plot_pumpdirection
    case 'Vertical' 
        xlim(plotaxis,[minWL,maxWL]);
    case 'Horizontal'
        ylim(plotaxis,[minWL,maxWL]);
    end
else
    switch plot_pumpdirection
    case 'Vertical' 
        ylim(plotaxis,[minWL,maxWL]);
    case 'Horizontal'
        xlim(plotaxis,[minWL,maxWL]);
    end
end


% --- Executes during object creation, after setting all properties.
function editWLmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function plot_Ncontours_Callback(hObject, eventdata, handles)
handles = UpdateEdits_2DIR(handles,0);
plot_2DIR(handles,handles.MainAxes);


% --- Executes during object creation, after setting all properties.
function plot_Ncontours_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editMaxZ_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function editMaxZ_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in timescalemenu.
function timescalemenu_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function timescalemenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_Nwhites_Callback(hObject, eventdata, handles)
handles = UpdateEdits_2DIR(handles,0);
plot_2DIR(handles,handles.MainAxes);


% --- Executes during object creation, after setting all properties.
function plot_Nwhites_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ShowContoursTick.
function ShowContoursTick_Callback(hObject, eventdata, handles)
handles = UpdateEdits_2DIR(handles,0);
plot_2DIR(handles,handles.MainAxes);



function plot_contourLineStyle_Callback(hObject, eventdata, handles)
handles = UpdateEdits_2DIR(handles,0);
plot_2DIR(handles,handles.MainAxes);


% --- Executes during object creation, after setting all properties.
function plot_contourLineStyle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plot_pumpdirection.
function plot_pumpdirection_Callback(hObject, eventdata, handles)
handles = UpdateEdits_2DIR(handles,0);
plot_2DIR(handles,handles.MainAxes);



% --- Executes during object creation, after setting all properties.
function plot_pumpdirection_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plot_axislegend.
function plot_axislegend_Callback(hObject, eventdata, handles)
handles = UpdateEdits_2DIR(handles,0);
plot_2DIR(handles,handles.MainAxes);



% --- Executes during object creation, after setting all properties.
function plot_axislegend_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot_colorrange_Callback(hObject, eventdata, handles)
handles = UpdateEdits_2DIR(handles,0);
plot_2DIR(handles,handles.MainAxes);

% --- Executes during object creation, after setting all properties.
function plot_colorrange_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot_skiplevels_Callback(hObject, eventdata, handles)
handles = UpdateEdits_2DIR(handles,0);
plot_2DIR(handles,handles.MainAxes);


% --- Executes during object creation, after setting all properties.
function plot_skiplevels_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Population_delay.
function Population_delay_Callback(hObject, eventdata, handles)
handles = UpdateEdits_2DIR(handles,0);
plot_2DIR(handles,handles.MainAxes);
if handles.DataTypeMenu.Value ~= 4
    td = handles.ShowTimeDomain.Value;
    ph = handles.ShowPhasing.Value;
    pc = handles.ShowProbeCalibration.Value;
    if td == 1
        what = 'td';
    elseif ph == 1
        what = 'ph';
    elseif pc == 1
        what = 'pc';
    end
    SecondaryPlot(handles,handles.SecondaryAxes,what);
end


% --- Executes during object creation, after setting all properties.
function Population_delay_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_contourfill.
function plot_contourfill_Callback(hObject, eventdata, handles)
handles = UpdateEdits_2DIR(handles,0);
plot_2DIR(handles,handles.MainAxes);


% --- Executes on selection change in plot_colourscheme.
function plot_colourscheme_Callback(hObject, eventdata, handles)
handles = UpdateEdits_2DIR(handles,0);
plot_2DIR(handles,handles.MainAxes);


% --- Executes during object creation, after setting all properties.
function plot_colourscheme_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in SaveProbeCal.
function SaveProbeCal_Callback(hObject, eventdata, handles)
ProbeAxis = handles.ProbeAxis;
rootdir = handles.rootdir;
probecalfile= 'CalibratedProbe.csv';
if exist(probecalfile,'file') == 0
    csvwrite([rootdir filesep probecalfile],ProbeAxis);
else
    % Ask the user what to do
    choice = questdlg('A calibrated probe file already exists, overwrite?', ...
        'Probe calibration', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            csvwrite([rootdir filesep probecalfile],ProbeAxis);
        case 'No'
            % Do nothing
    end
end

% --- Executes on button press in SymColRange_tick.
function SymColRange_tick_Callback(hObject, eventdata, handles)
handles = UpdateEdits_2DIR(handles,0);
plot_2DIR(handles,handles.MainAxes);

% --- Executes on button press in CutPlot_tick.
function CutPlot_tick_Callback(hObject, eventdata, handles)
handles = UpdateEdits_2DIR(handles,0);
plot_2DIR(handles,handles.MainAxes);


function phase_coeffs_Callback(hObject, eventdata, handles)
% This shouldn't do anything

function phase_coeffs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GoTodayDir.
function GoTodayDir_Callback(hObject, eventdata, handles)
% Go to the folder
datedir     = datestr(now,'yyyymmdd');
todaydir    = [handles.defaultdir filesep datedir];
if exist(todaydir,'dir') ~= 0
    cd (todaydir);
    % Update the string
    handles.CurrDir.String = pwd;
    handles.rootdir = handles.CurrDir.String;
    % Put list of subfolders into subFolders
        % Get a list of all files and folders in this folder.
        files = dir(handles.rootdir);
        % Extract only those that are directories.
        subFolders = files([files.isdir]);
    % Go to the selected folder
        cd (handles.rootdir);
    % Sort the names and create variables for the list
        dir_struct = subFolders;
        [sorted_names,sorted_index] = sortrows({dir_struct.name}');
        handles.file_names = sorted_names;
        handles.is_dir = [dir_struct.isdir];
        handles.sorted_index = sorted_index;
        guidata(hObject,handles)
    % Update DatafoldersList, removing ".." and "." from the visible list
        handles.file_names = handles.file_names(3:end);
        set(handles.DatafoldersList,'String',handles.file_names,'Value',1)
    guidata(hObject,handles)
end

function EditProbeAxis_tick_Callback(hObject, eventdata, handles)
% This shouldn't do anything

% --- Executes on selection change in Transient2D_mode.
function Transient2D_mode_Callback(hObject, eventdata, handles)
% Reprocess the data without reloading
handles = process2DIR(handles,1);
guidata(hObject,handles);
% Clear error string
set(handles.ErrorText,'String', "");
% Update the edit controls
handles = UpdateEdits_2DIR(handles,1);
% Plot the data in the main window (for preview)
plot_2DIR(handles,handles.MainAxes);   
set(handles.MainAxes,'Visible','On')
% Generate the secondary plot (show phasing data by default)
SecondaryPlot(handles,handles.SecondaryAxes,'ph');
set(handles.SecondaryAxes,'Visible','On')
    handles.ShowTimeDomain.Value = 0;
    handles.ShowProbeCalibration.Value = 0;
    handles.ShowPhasing.Value = 1;
% Delete progress bar
delete(handles.WaitBar);


% --- Executes during object creation, after setting all properties.
function Transient2D_mode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ShowSpecDiff_Callback(hObject, eventdata, handles)
handles = UpdateEdits_2DIR(handles,0);
plot_2DIR(handles,handles.MainAxes);

% --- Executes on button press in SubtractSpectra.
function SubtractSpectra_Callback(hObject, eventdata, handles)
% Ask the user what to plot
    sel_spectra     = inputdlg({'Enter the index of dataset A:','Enter the index of dataset B:'},'Subtract spectra (A-B)',[1 50],{'1','2'});
    if isempty(sel_spectra)
        return
    end

% Get the indices of the selected spectra
    A_idx = str2double(sel_spectra{1});
    B_idx = str2double(sel_spectra{2});

% Process the spectra
    % Create a copy of the current handles structure to load the second dataset
    handles_B   = handles;
    file_list   = handles.DatafoldersList.String;
    
    %%% Process A
        % Define filename
        handles.datafilename  = string(cellstr(file_list{A_idx})); 
        % Load data
        handles               = load2DIRlab1(handles);
        % Phase spectra
        handles               = process2DIR(handles,0);
        % Delete progress bar
        delete(handles.WaitBar);

    %%% Process B
        % Define filename
        handles_B.datafilename  = string(cellstr(file_list{B_idx})); 
        % Load data
        handles_B               = load2DIRlab1(handles_B);
        % Phase spectra
        handles_B               = process2DIR(handles_B,0);
        % Delete progress bar
        delete(handles_B.WaitBar);
        
% Calculate the differential 2D-IR spectra
    % Get the phased 2D data
    A_data  = handles.PROC_2D_DATA;
    B_data  = handles_B.PROC_2D_DATA;
    
    % Check if the datasets have the number of t2 delays
    if length(A_data) ~= length(B_data)
        opts.Interpreter    = 'tex';
        opts.WindowStyle    = 'replace';
        warndlg('Datasets do not have the same number of \itt_2\rm delays!','Warning',opts);
        return
    end
  
    % Check if the datasets have the same FT and pixel size
    if length(A_data{1,1}) ~= length(B_data{1,1})
        opts.Interpreter    = 'tex';
        opts.WindowStyle    = 'replace';
        warndlg('Datasets do not have the same \omega_1 size!','Warning',opts);
        return
    end
    
    % Preallocate variables
    Ncols       = size(A_data,2);
    diff_data   = cell(length(A_data),Ncols);
    
    % If everything is OK, do the subtraction
    for i=1:length(A_data)
        diff_data{i}    = A_data{i} - B_data{i};
    end
    
% Write back to main handles structure
    handles.PROC_2D_DATA = diff_data;

% Update controls, plot the data and save to the main handles structure
    % Update the edit controls
    handles = UpdateEdits_2DIR(handles,1);
    % Plot the data in the main window (for preview)
    plot_2DIR(handles,handles.MainAxes);   
    set(handles.MainAxes,'Visible','On')
    % Generate the secondary plot (show phasing data by default)
    SecondaryPlot(handles,handles.SecondaryAxes,'ph');
    set(handles.SecondaryAxes,'Visible','On')
        handles.ShowTimeDomain.Value = 0;
        handles.ShowProbeCalibration.Value = 0;
        handles.ShowPhasing.Value = 1;
    % Enable disabled controls
    EnableControls_2DIR(handles)
    % Update handles
    guidata(hObject,handles)
