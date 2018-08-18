function varargout = DataAnalysis_GUI(varargin)
% DATAANALYSIS_GUI MATLAB code for DataAnalysis_GUI.fig
%      DATAANALYSIS_GUI, by itself, creates a new DATAANALYSIS_GUI or raises the existing
%      singleton*.
%
%      H = DATAANALYSIS_GUI returns the handle to a new DATAANALYSIS_GUI or the handle to
%      the existing singleton*.
%
%      DATAANALYSIS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATAANALYSIS_GUI.M with the given input arguments.
%
%      DATAANALYSIS_GUI('Property','Value',...) creates a new DATAANALYSIS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DataAnalysis_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DataAnalysis_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DataAnalysis_GUI

% Last Modified by GUIDE v2.5 07-Aug-2018 21:45:23

% Ricardo Fernández-Terán - v3.9a - 11.06.2018
% ----CHANGELOG:
% * Fixed a bug when plotting normalised kinetics per scan (wasn't dividing by the abs value)
% * Updated all plotting routines to make plot format consistent. Font sizes and styles should now be the same in all plots.
% * Included loading routine for TRES data files (*.dat)
% * Updated tWLshift so it doesn't overwrite the .bak file, only _wavenumbers.csv is overwritten!
% * Included loading routines for Lab 4 data (old format)
% 
% ----PENDING:
% * Add option to subtract solvent pFID and/or scattering surface (array - array) to get the "clean" signals

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DataAnalysis_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @DataAnalysis_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

% Disable annoying warnings
warning('off','MATLAB:Axes:NegativeLimitsInLogAxis');
warning('off','MATLAB:legend:IgnoringExtraEntries');

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before DataAnalysis_GUI is made visible.
function DataAnalysis_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure,
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DataAnalysis_GUI (see VARARGIN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default directory:
if exist('C:\GUIoptions.txt','file')==2
    handles.defaultdir = readParam('defaultdir','C:\GUIoptions.txt');
else
    handles.defaultdir = [];
end
% Make static text renderable using LaTex
    % TEXT annotations need an axes as parent so create an invisible axes which
    % is as big as the figure
    handles.laxis = axes('parent',hObject,'units','normalized','position',[0 0 1 1],'visible','off');
    % Find all static text UICONTROLS whose 'Tag' starts with latex_
    lbls = findobj(hObject,'-regexp','tag','latex_*');
    for i=1:length(lbls)
          l = lbls(i);
          % Get current text, position and tag
          set(l,'units','normalized');
          s = get(l,'string');
          p = get(l,'position');
          t = get(l,'tag');
          % Remove the UICONTROL
          delete(l);
          % Replace it with a TEXT object 
          handles.(t) = text(p(1),p(2),s,'interpreter','latex');
    end
% Choose default command line output for DataAnalysis_GUI
handles.output = hObject;
% Linear time axis by default
handles.linlog = 'lin';
handles.linlogtick.Value=0;
% Symmetric colour range by default
handles.SymmetricColourRange_tick.Value=1;
% At the beginning there is no FTIR spectrum loaded
handles.IRloaded = 0;
% By default, do not save kinetic traces
handles.DoSaveTraces = 0;

% UIWAIT makes DataAnalysis_GUI wait for user response (see UIRESUME)
% uiwait(handles.MainGUI); DO NOT USE!

% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = DataAnalysis_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
%varargout{1} = handles.output;


% --- Executes on button press in BrowseRootDirButton.
function BrowseRootDirButton_Callback(hObject, eventdata, handles)
% Ask user to select root directory (to load all data subdirectories
% in this directory)
rootdir = uigetdir(handles.defaultdir);
handles.rootdir = rootdir;

switch handles.DataTypeMenu.Value
    case 4 % *.dat files
        if rootdir ~= 0
            % Update text accordingly
            handles.CurrDir.String  = rootdir;
            % Get a list of all files and folders in this folder.
            files = dir(rootdir);
            % Get only those that are *.dat
            files = files(contains({files.name},".dat",'IgnoreCase',1));
            % Sort the names and create variables for the list
                [sorted_names,sorted_index] = sortrows({files.name}');
                handles.file_names          = sorted_names;
                handles.sorted_index        = sorted_index;
                guidata(hObject,handles)
                set(handles.DatafoldersList,'String',handles.file_names,'Value',1)
                handles.is_dir = 0;
                guidata(hObject,handles)
        else 
            handles.CurrDir.String = "Please select a directory!  :(";
        end
    otherwise % = 1,2,3
        if rootdir ~= 0
            % Update text accordingly
            handles.CurrDir.String  = rootdir;
            % Get a list of all files and folders in this folder.
            files = dir(rootdir);
            % Extract only those that are directories.
            subFolders = files([files.isdir]);
%           % Go to the folder
%           cd (rootdir);
            % Sort the names and create variables for the list
                [sorted_names,sorted_index] = sortrows({subFolders.name}');
                handles.file_names = sorted_names;
                handles.is_dir = [subFolders.isdir];
                handles.sorted_index = sorted_index;
                guidata(hObject,handles)
            % Update DatafoldersList, removing ".." and "." from the visible list
                handles.file_names(1:2,:) = [];
                set(handles.DatafoldersList,'String',handles.file_names,'Value',1)
                guidata(hObject,handles)     
        %       else
        %         % LOAD DATA IN LAB 4 PUMP-PROBE FORMAT
        %         % Execute PatchFileList_Lab4 and use the output to build the file list
        %             handles = PatchFileList_Lab4(handles,'1D');
        %             guidata(hObject,handles);
        %             set(handles.DatafoldersList,'String',handles.SampleName,'Value',1);
        %       end
        else 
            set(handles.CurrDir, 'String', 'Please select a directory!  :(');
        end
end

% --- Executes on selection change in DatafoldersList.
function DatafoldersList_Callback(hObject, eventdata, handles)
switch get(gcf,'SelectionType')
case 'normal' % Single click
    index_selected  = handles.DatafoldersList.Value;
    file_list       = handles.DatafoldersList.String;
    datafilename    = cellstr(file_list{index_selected}); % Item selected in list box
    
    % Write datafilename to handles and update
    handles.datafilename = string(datafilename);
    guidata(hObject,handles)
    
% If directory
    if  handles.is_dir(handles.sorted_index(index_selected))
        % Load the data (replace handles) and update the main handles
        error=0; handles.Nscans = NaN;
        try
           switch handles.DataTypeMenu.Value
               case 1 % Lab 2: UV-Vis Pump - IR Probe (ns)
                   handles = LoadDataIR2(handles);
                   tempdir = [char(handles.CurrDir.String) filesep char(datafilename) filesep 'temp'];
                   if exist(tempdir,'dir') == 7
                       cd(tempdir);
                       filelist = dir('*signal*.csv');
                       handles.Nscans = floor(length(filelist)./2);
                       handles.Nscans_number.String = num2str(handles.Nscans);
                   else
                       handles.Nscans_number.String = 'N/A';
                       handles.KineticsPerScan.Enable      = 'Off';
                       handles.SpectraPerScan.Enable       = 'Off';
                       handles.BinScans.Visible            = 'Off';
                       handles.BinScans.Enable             = 'Off';
                       handles.text35.Visible              = 'Off';
                       handles.text36.Visible              = 'Off';
                       handles.text35.Enable               = 'Off';
                       handles.text36.Enable               = 'Off';
                   end
                   set(handles.ErrorText,'String', '');
                   handles.SymmetricColourRange_tick.Value = 1;
               case 2 % Lab 1: UV-Vis/IR Pump - IR probe (fs)
                   handles  = LoadDataIR1(handles,1);
                   handles.SlowMod_selector.Value   = 1;
                   % TEMP dir story
                   tempdir         = [char(handles.CurrDir.String) filesep char(datafilename) filesep 'temp'];
                   tempdirOUT      = [char(handles.CurrDir.String) filesep char(datafilename) 'temp'];
                   % Check if the temp dir is outside (Lab 1) or inside the data folder (Lab 4)
                   if exist(tempdir,'dir') == 0
                       tempdir     = tempdirOUT;
                   end
                   % Check if the temp dir really exists
                   if exist(tempdir,'dir') == 7
                       handles.Nscans_number.String = num2str(handles.Nscans);
                   else
                       handles.Nscans_number.String        = 'N/A';
                       handles.KineticsPerScan.Enable      = 'Off';
                       handles.SpectraPerScan.Enable       = 'Off';
                       handles.BinScans.Visible            = 'Off';
                       handles.BinScans.Enable             = 'Off';
                       handles.text35.Visible              = 'Off';
                       handles.text36.Visible              = 'Off';
                       handles.text35.Enable               = 'Off';
                       handles.text36.Enable               = 'Off';
                   end
                   % Clear error string
                   set(handles.ErrorText,'String', '');
                   handles.SymmetricColourRange_tick.Value = 1;
           end
        catch err
           set(handles.ErrorText,'String', string(err.message));
           error=1;
        end
        if error == 0
                set(handles.axes1,'Visible','On')
                % Initialise defaults and update the Edits for the first time:
                handles = UpdateEdits(handles);
                guidata(hObject,handles)
                % Override default w/current status of linlog tick
                switch handles.linlogtick.Value
                    case 0
                        handles.linlog = 'lin';
                    case 1
                        handles.linlog = 'log';
                end
        % Plot the data in the main window (for preview)
        % Takes into account the current status of BkgSubTick
        if exist('handles.rawcorr','var') == 0
            switch handles.BkgSubTick.Value
                case 1
                    handles.rawcorr='CORRECTED';
                case 0
                    handles.rawcorr='RAW';
            end
        end
        switch handles.rawcorr
            case 'RAW'
                plot2D(handles,handles.rawsignal,handles.axes1,'On');
            case 'CORRECTED'
                plot2D(handles,handles.corrdata,handles.axes1,'On');
        end
        % Update handles
        guidata(hObject,handles)
        end
% If NOT a dir        
    else 
        % Load the data (replace handles) and update the main handles
        error=0; handles.Nscans = NaN;
        try
           switch handles.DataTypeMenu.Value
               case 4 % TRES (*.dat)
                   handles = LoadDataTRES(handles);
                   handles.SymmetricColourRange_tick.Value = 0;
               case 3 % Lab 4 data (old format)
                   handles = LoadDataIR4(handles);
           end
        catch err
           set(handles.ErrorText,'String', string(err.message));
           error=1;
        end
        if error == 0
                set(handles.axes1,'Visible','On')
                % Initialise defaults and update the Edits for the first time:
                handles = UpdateEdits(handles);
                guidata(hObject,handles)
                % Override default w/current status of linlog tick
                switch handles.linlogtick.Value
                    case 0
                        handles.linlog = 'lin';
                    case 1
                        handles.linlog = 'log';
                end
        % Plot the data in the main window (for preview)
        % Takes into account the current status of BkgSubTick
        if exist('handles.rawcorr','var') == 0
            switch handles.BkgSubTick.Value
                case 1
                    handles.rawcorr='CORRECTED';
                case 0
                    handles.rawcorr='RAW';
            end
        end
        switch handles.rawcorr
            case 'RAW'
                plot2D(handles,handles.rawsignal,handles.axes1,'On');
            case 'CORRECTED'
                plot2D(handles,handles.corrdata,handles.axes1,'On');
        end
        % Update handles
        guidata(hObject,handles)
        end
    end
case 'open'
        % Go to the selected folder
        index_selected  = handles.DatafoldersList.Value;
        
    if handles.is_dir(handles.sorted_index(index_selected))
        file_list       = handles.DatafoldersList.String;
        datafilename    = cellstr(file_list{index_selected}); % Item selected in list box
        % Write datafilename to handles(datafilename) and update
        dirnamech       = char(datafilename);
        % Go to the selected dir
        currdir         = char(handles.CurrDir.String);
        cd([currdir filesep dirnamech]);
        % Update the string
        handles.CurrDir.String      = pwd;
        handles.rootdir             = handles.CurrDir.String;
        % Put list of subfolders into subFolders
            % Get a list of all files and folders in this folder.
            files = dir(handles.rootdir);
            % Get a logical vector that tells which is a directory.
            dirFlags = [files.isdir];
            % Extract only those that are directories.
            subFolders = files(dirFlags);
        % Sort the names and create variables for the list
            dir_struct = subFolders;
            [sorted_names,sorted_index] = sortrows({dir_struct.name}');
            handles.file_names = sorted_names;
            handles.is_dir = [dir_struct.isdir];
            handles.sorted_index = sorted_index;
            guidata(hObject,handles)
        % Update DatafoldersList, removing ".." and "." from the visible list
            handles.file_names(1:2,:) = [];
            set(handles.DatafoldersList,'String',handles.file_names,'Value',1)
            guidata(hObject,handles)
    end
end


% --- Executes during object creation, after setting all properties.
function DatafoldersList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DatafoldersList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlotBkg.
function PlotBkg_Callback(hObject, eventdata, handles)
% Create a new figure with consistent format
fh = figure();
fh.Position(3)  = 800;
fh.Position(4)  = 425;
fh.Color        = [1 1 1];
% Define the axes
handles.axes2 = axes('Parent',fh);
axes(handles.axes2);
handles.axes2.Units     = 'pixels';
handles.axes2.Position  = [80 65 680 320];
% Do the plot
plot(handles.cmprobe,handles.bkg,'LineWidth',1.5,'Color','r')
% Formatting
axis tight
set(gca,'FontSize',12);
xlabel(['Wavenumbers',' (cm^{-1})'],'FontSize',13,'FontWeight','bold');
ylabel('Background signal (m\DeltaOD)','FontSize',13,'FontWeight','bold');
title(['Background of ',char(handles.datafilename)],'Interpreter','none','FontSize',12);
hline = refline(0,0); hline.Color = [0.5 0.5 0.5];

% --- Executes on button press in make2Dplot.
function make2Dplot_Callback(hObject, eventdata, handles)
linlog  = handles.linlog;
fh      = figure();
newaxes = axes('Parent',fh);
switch handles.rawcorr
    case 'RAW'
        plot2D(handles,handles.rawsignal,newaxes,'Off')
    case 'CORRECTED'
        plot2D(handles,handles.corrdata,newaxes,'Off')
end
set(newaxes,'yscale',linlog)

% Make the plot format constant
fh.Units='pixels';
fh.Position(2)=fh.Position(2).*0.7;
fh.Position(3)=590; % Width in pixels
fh.Position(4)=510; % Height in pixels
fh.Color = [1 1 1];
newaxes.FontSize = 12;
newaxes.Units = 'pixels';
newaxes.Position = [75 75 420 385];
% Make uniform, consistent format
newaxes.LineWidth = 1;
newaxes.TickLength = [0.015 0.035];

% Make the axes resizable
newaxes.Units = 'normalized';

guidata(hObject,handles)

% --- Executes on button press in make3Dplot.
function make3Dplot_Callback(hObject, eventdata, handles)
% hObject    handle to make3Dplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
linlog = handles.linlog;
fh = figure();
fh.Color = [1 1 1];
newaxes = axes('Parent',fh);
switch handles.rawcorr
    case 'RAW'
        plot3D(handles,handles.rawsignal,newaxes)
    case 'CORRECTED'
        plot3D(handles,handles.corrdata,newaxes)
end
set(handles.newaxes,'yscale',linlog)

% Make the plot format constant
fh.Units='pixels';
fh.Position(2)=fh.Position(2).*0.7;
fh.Position(3)=590; % Width in pixels
fh.Position(4)=510; % Height in pixels
fh.Color = [1 1 1];
newaxes.FontSize = 12;
newaxes.Units = 'pixels';
newaxes.Position = [75 75 420 385];
% Make uniform, consistent format
newaxes.LineWidth = 1;
newaxes.TickLength = [0.015 0.035];

% Make the axes resizable
newaxes.Units = 'normalized';

guidata(hObject,handles)


% --- Executes on button press in doSVD.
function doSVD_Callback(hObject, eventdata, handles)
% hObject    handle to doSVD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch handles.rawcorr
    case 'RAW'
        data = handles.rawsignal;
    case 'CORRECTED'
        data = handles.corrdata;
end
[T,SV,WL] = svd(data);
handles.SVD = 1;
n = str2num(handles.nSVD.String);
% Plot the Singular Values
figure()
sv_axes = subplot(2,2,[3 4]);
cmap = colormap(lines(length(n)));
SVAL = diag(SV); SVAL = SVAL(n,:); SV_diag=diag(SVAL,0);
scatter(n,SVAL,75,cmap,'filled','LineWidth',0.75);
axis tight; yl = ylim; ylim([yl(1) 1.1.*yl(2)]); set(gca, 'FontSize', 10); set(gca,'xtick',n)
%title(['Singular values of ',handles.datafilename],'Interpreter','none','FontWeight','bold');
xlabel('SVD index','FontWeight','bold');
ylabel('SVD magnitude','FontWeight','bold');
% Plot the RH vectors (Wavelength)
subplot(2,2,1)
plot(handles.cmprobe,WL(:,n)); axis tight;
xlabel('Wavenumbers (cm^{-1})','FontWeight','bold');
ylabel('\DeltaAbs (mOD)','FontWeight','bold');
% Plot the LH vectors (Time)
subplot(2,2,2)
plot(handles.delays,T(:,n)*SV_diag); axis tight;
xlabel(['Delays ( ',handles.timescale,')'],'FontWeight','bold');

%Add the cumulative information content as function of increasing SVD
%components
axes(sv_axes)
yyaxis(sv_axes,'right')
sv_axes.YColor = 'r';
cumulSV=zeros(length(n),1);
for q=1:length(n)
    cumulSV(q,1)=sum(SVAL(1:q))/sum(SVAL)*100;
end
plot(transpose(n),cumulSV,'-or','LineWidth',0.75);
text(transpose(n),cumulSV,num2str(cumulSV,3),'HorizontalAlignment','left','VerticalAlignment','top')
ylabel('Cumulative content (%)','FontWeight','bold');
ylim([70 100])

% --- Executes on button press in SVD_fit.
function SVD_fit_Callback(hObject, eventdata, handles)
handles.FitType = 2;
switch handles.InteractiveModeTick.Value
    case 1
        % Read SVD components to fit from the control in the main window
        SingVals    = transpose(str2num(handles.nSVD.String));  
    case 0
        % Ask the user for the SVD components to fit
        prompt      = 'Enter desired SVD components to fit (i.e. 1 to 3 = 1:3):';
        dlg_title   = 'Select SVD components';
        num_lines   = [1 60];
        SelSVD      = inputdlg(prompt,dlg_title,num_lines);
        SingVals    = transpose(str2double(SelSVD{1}));
end
nSVD = length(SingVals);
% Determine whether the data is RAW or CORRECTED
switch handles.rawcorr
    case 'RAW'
        data=handles.rawsignal;
    case 'CORRECTED'
        data=handles.corrdata;
end
% Do the SVD and take only the components selected by the user
[T,SV,WL] = svd(data);
handles.SVD = 1;
SV_diag=diag(diag(SV),0);
% Prepare the labels and required variables
YLabel = '\DeltaAbs (mOD)';
FitType=2; % i.e. 1= LOCAL, 2 = GLOBAL SVD fit
% Now open the fit window
[fitP,fitP_SD,FitFunc,Residuals,FitDone,tau_index,c_index,InputParams] = FitWindow(handles.delays,T(:,SingVals)*SV_diag(SingVals,SingVals),SingVals,YLabel,handles.timescale,FitType);
% After the user has closed the window, save the global fit results to file
% (overwrite by default)
    %%%% TO BE DONE
% Get the Taus (and SD's)
    Taus            = fitP(tau_index);
    Taus_SD         = fitP_SD(tau_index);
    NumExp          = length(Taus);
    C_matrix        = fitP(c_index(:,1:NumExp));
    C_SD_matrix     = fitP_SD(c_index(:,1:length(Taus)));
% Plot the DAS
    % Get the DAS
    SelWL = WL(:,SingVals);
    DAS = SelWL*C_matrix;  
    % Create a new figure with consistent format
    fh = figure();
    fh.Position(3)  = 800;
    fh.Position(4)  = 425;
    fh.Color        = [1 1 1];

    % Define the axes
    handles.axes2 			= axes('Parent',fh);
    handles.axes2.Units     = 'pixels';
    handles.axes2.Position  = [75 75 680 320];
    axes(handles.axes2);
    % Plot the DAS
    cm=colormap(othercolor('Mrainbow',NumExp+1));
    legendlabel={};
    for n=1:NumExp
       % Check if Tau > 1000, then convert to next prefix
       if Taus(n) > 1000
           switch handles.timescale
                case 'ns'
                    newtimescale = ['\mu' 's'];
                case 'ps'
                    newtimescale = 'ns';
                case 'fs'
                    newtimescale = 'ps';
           end
           f=1000;
       elseif Taus(n) < 1
           switch handles.timescale
                case 'ns'
                    newtimescale = 'ps';
                case 'ps'
                    newtimescale = 'fs';
           end
           f=1/1000;
       else
           f=1;
           newtimescale = handles.timescale;
       end
       strTau      = num2str(Taus(n)/f,'%.3g');
       strTau_SD   = num2str(Taus_SD(n)/f,'%.3g');
       
       legendlabel{n}=['\tau_{' num2str(n) '}' ' = ' strTau ' \pm ' strTau_SD ' ' newtimescale];
       plot(handles.cmprobe,DAS(:,n),'LineWidth',2,'Marker','none','color',cm(n,:));
       hold on
    end
    % Calculate the infinite component (if it exists)
    if prod(InputParams.C_holds(:,5)) == 0
        cInfty          = fitP(c_index(:,5));
        cInfty_SD       = fitP_SD(c_index(:,5));
        INF = SelWL*cInfty;
        % Plot the infinite component
        plot(handles.cmprobe,INF,'LineWidth',2,'Marker','none','color',cm(end,:));
        legendlabel{NumExp+1}='\tau_{\infty}';
    elseif sum(InputParams.C_holds(:,5)) ~= 0
        cInfty          = InputParams.C_values(:,6);
    end
    % Calculate the offset (if it exists)
    if sum(InputParams.C_holds(:,6)) ~= nSVD && sum(InputParams.C_values(:,6)) ~= 0 
        Offset          = fitP(c_index(:,6));
        Offset_SD       = fitP_SD(c_index(:,6));
    elseif sum(InputParams.C_values(:,6)) ~= 0
        Offset          = InputParams.C_values(:,6);
    end
    
    % Formatting
    l=legend(legendlabel); legend('boxoff'); legend('Location','northeast');
    l.Interpreter = 'tex';
    set(gca,'FontSize',12);
    xlabel('Wavenumbers (cm^{-1})','FontSize',13,'FontWeight','bold')
    ylabel('Amplitude','FontSize',13,'FontWeight','bold')
    axis tight
    hline = refline(0,0); hline.Color = [0.5 0.5 0.5];
    set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

% Plot the 2D residual surface
    % Initialise variables
    FitEval=zeros(length(handles.delays),nSVD);
    GlobalFit=zeros(length(handles.delays),length(handles.cmprobe));
    GlobalRes=zeros(length(handles.delays),length(handles.cmprobe));
    % Evaluate the fit functions with the fitted params
    for m=1:nSVD
        FitEval(:,m) = feval(FitFunc{m},fitP,handles.delays);
    end
    % Retrieve the 2D maps (and residuals)
    GlobalFit=FitEval*transpose(SelWL);
    GlobalRes=data-GlobalFit;
    % Plot the fitted surface
    figure;
    ax = axes();
    plot2D(handles,GlobalFit,ax,'Off');
    title({handles.datafilename;['FITTED SURFACE';'']},'Interpreter','none')
    % Plot the residual surface
    figure
    ax = axes();
    plot2D(handles,GlobalRes,ax,'Off')
    title({handles.datafilename;['RESIDUAL SURFACE';'']},'Interpreter','none');
    % Plot the residual surface with maximum scale
    figure
    ax = axes();
    plot2D(handles,GlobalRes,ax,'Off')
    caxis auto
    title({handles.datafilename;['RESIDUAL SURFACE - max Z scale';'']},'Interpreter','none');   
% Everything Done! :)

function editWLmax_Callback(hObject, eventdata, handles)
% hObject    handle to editWLmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plotranges(4) = str2double(get(hObject,'String'));
guidata(hObject,handles)
WLlim = handles.plotranges(3:4);
axes(handles.axes1)
try
    xlim(WLlim)
catch err
    errordlg(err.message,'Expression Error');
end
guidata(hObject,handles)
% Hints: get(hObject,'String') returns contents of editWLmax as text
%        str2double(get(hObject,'String')) returns contents of editWLmax as a double


% --- Executes during object creation, after setting all properties.
function editWLmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWLmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editWLmin_Callback(hObject, eventdata, handles)
% hObject    handle to editWLmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plotranges(3) = str2double(get(hObject,'String'));
guidata(hObject,handles)
WLlim = handles.plotranges(3:4);
axes(handles.axes1)
try
    xlim(WLlim)
catch err
    errordlg(err.message,'Expression Error');
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function editWLmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editTmax_Callback(hObject, eventdata, handles)
handles.plotranges(2) = str2double(get(hObject,'String'));
guidata(hObject,handles)
delaylim = handles.plotranges(1:2);
axes(handles.axes1)
try
    ylim(delaylim)
catch err
    errordlg(err.message,'Expression Error');
end
guidata(hObject,handles)
% Hints: get(hObject,'String') returns contents of editTmax as text
%        str2double(get(hObject,'String')) returns contents of editTmax as a double


% --- Executes during object creation, after setting all properties.
function editTmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editTmin_Callback(hObject, eventdata, handles)
% hObject    handle to editTmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plotranges(1) = str2double(get(hObject,'String'));
guidata(hObject,handles)
delaylim = handles.plotranges(1:2);
axes(handles.axes1)
try
    ylim(delaylim);
catch err
    errordlg(err.message,'Expression Error');
end
guidata(hObject,handles)
% Hints: get(hObject,'String') returns contents of editTmin as text
%        str2double(get(hObject,'String')) returns contents of editTmin as a double


% --- Executes during object creation, after setting all properties.
function editTmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNcontours_Callback(hObject, eventdata, handles)
% hObject    handle to editNcontours (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% TO FIX: MAKE IT CHANGE THE PLOT WITHOUT DOING IT AGAIN
handles.plotranges(6) = str2double(get(hObject,'String'));
switch handles.rawcorr
    case 'RAW'
        plot2D(handles,handles.rawsignal,handles.axes1,'On')
    case 'CORRECTED'
        plot2D(handles,handles.corrdata,handles.axes1,'On')
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function editNcontours_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNcontours (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMaxZ_Callback(hObject, eventdata, handles)
% % Zlimits = str2num(get(hObject,'String'));
% % axes(handles.axes1)
% % if size(Zlimits,2) == 1
% %     handles.plotranges(5) = Zlimits;
% %     Zminmax = handles.plotranges(5);
% %     caxis([-Zminmax,Zminmax])
% % elseif size(Zlimits,2) == 2
% %     handles.plotranges(5) = max(abs(Zlimits));
% %     caxis(Zlimits)
% %     
% % end
switch handles.rawcorr
    case 'RAW'
        plot2D(handles,handles.rawsignal,handles.axes1,'On')
    case 'CORRECTED'
        plot2D(handles,handles.corrdata,handles.axes1,'On')
end
guidata(hObject,handles)




% --- Executes during object creation, after setting all properties.
function editMaxZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)


% --- Executes on selection change in timescalemenu.
function timescalemenu_Callback(hObject, eventdata, handles)
% hObject    handle to timescalemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get the selection
timescale = cellstr(get(hObject,'String'));
% Update handles and text on GUI
handles.timescale = timescale{get(hObject,'Value')};
set(handles.unitsBkgsub,'String',handles.timescale);
set(handles.tmin_text,'String',strcat('Min time (',handles.timescale,')'));
set(handles.tmax_text,'String',strcat('Max time (',handles.timescale,')'));
axes(handles.axes1);
ylabel(['Delays',' (',handles.timescale,')'])
guidata(hObject,handles)
% Hints: contents = cellstr(get(hObject,'String')) returns timescalemenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from timescalemenu


% --- Executes during object creation, after setting all properties.
function timescalemenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timescalemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in linlogtick.
function linlogtick_Callback(hObject, eventdata, handles)
% hObject    handle to linlogtick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(hObject,'Value')
    case 1
        set(handles.axes1,'yscale','log')
        handles.linlog = 'log';
    case 0
        set(handles.axes1,'yscale','lin')
        handles.linlog = 'lin';
end
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of linlogtick


% --------------------------------------------------------------------
function Zoom_In_OnCallback(hObject, eventdata, handles)
% hObject    handle to Zoom_In (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plotranges(1:2) = ylim;
handles.plotranges(3:4) = xlim;
guidata(hObject,handles);
handles = UpdateEdits(handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function Zoom_In_OffCallback(hObject, eventdata, handles)
% hObject    handle to Zoom_In (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plotranges(1:2) = ylim;
handles.plotranges(3:4) = xlim;
guidata(hObject,handles);
handles = UpdateEdits(handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function Zoom_Out_OffCallback(hObject, eventdata, handles)
% hObject    handle to Zoom_Out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plotranges(1:2) = ylim;
handles.plotranges(3:4) = xlim;
guidata(hObject,handles);
handles = UpdateEdits(handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function Zoom_Out_OnCallback(hObject, eventdata, handles)
% hObject    handle to Zoom_Out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plotranges(1:2) = ylim;
handles.plotranges(3:4) = xlim;
handles = UpdateEdits(handles);
guidata(hObject,handles);



function mintimeBkg_Callback(hObject, eventdata, handles)
% hObject    handle to mintimeBkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mintimeBkg = str2double(get(hObject,'String'));
j = length(handles.delays);
t = handles.delays(j);
if mintimeBkg >= min(handles.delays)
    while t > mintimeBkg && j > 1
        j = j-1;
        t = handles.delays(j-1);
    end
    j = j-1;
    handles.j = j;
    k = handles.k;
    set(handles.mintimeBkg,'String',num2str(handles.delays(j)));
    set(handles.maxtimeBkg,'String',num2str(handles.delays(k)));
    guidata(hObject,handles);
    % Do the background subtraction
    handles.bkg = mean(handles.rawsignal(j:k,:));
    handles.corrdata = handles.rawsignal - handles.bkg;
    % Replot
    switch handles.rawcorr
        case 'RAW'
            plot2D(handles,handles.rawsignal,handles.axes1,'On')
        case 'CORRECTED'
            plot2D(handles,handles.corrdata,handles.axes1,'On')
    end
    guidata(hObject,handles);
else
    j = 1;
    handles.j = j;
    k = handles.k;
    set(handles.mintimeBkg,'String',num2str(handles.delays(j)));
    set(handles.maxtimeBkg,'String',num2str(handles.delays(k)));
    % Do the background subtraction
    handles.bkg = mean(handles.rawsignal(j:k,:));
    handles.corrdata = handles.rawsignal - handles.bkg;
    % Replot
    switch handles.rawcorr
        case 'RAW'
            plot2D(handles,handles.rawsignal,handles.axes1,'On')
        case 'CORRECTED'
            plot2D(handles,handles.corrdata,handles.axes1,'On')
    end
    guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function mintimeBkg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mintimeBkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxtimeBkg_Callback(hObject, eventdata, handles)
% hObject    handle to maxtimeBkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% TO FIX =  CHECK THAT THE TIME LIMITS ARE READ AND EXECUTED CORRECTLY
maxtimeBkg = str2double(get(hObject,'String'));
if maxtimeBkg <= max(handles.delays)
    k = 0;
    t = handles.delays(k+1);
    while t <= maxtimeBkg && k < (length(handles.delays)-1)
        k = k+1;
        t = handles.delays(k+1);
    end
    handles.k = k;
    j = handles.j;
    set(handles.mintimeBkg,'String',num2str(handles.delays(j)));
    set(handles.maxtimeBkg,'String',num2str(handles.delays(k)));
    guidata(hObject,handles);
    % Do the background subtraction and change status in handles.rawcorr
    handles.bkg = mean(handles.rawsignal(j:k,:));
    handles.corrdata = handles.rawsignal - handles.bkg;
    guidata(hObject,handles);
    % Replot
    switch handles.rawcorr
        case 'RAW'
            plot2D(handles,handles.rawsignal,handles.axes1,'On')
        case 'CORRECTED'
            plot2D(handles,handles.corrdata,handles.axes1,'On')
    end
    guidata(hObject,handles);
else
    k = length(handles.delays);
    handles.k = k;
    j = handles.j;
    set(handles.mintimeBkg,'String',num2str(handles.delays(j)));
    set(handles.maxtimeBkg,'String',num2str(handles.delays(k)));
    guidata(hObject,handles);
    % Do the background subtraction and change status in handles.rawcorr
    handles.bkg = mean(handles.rawsignal(j:k,:));
    handles.corrdata = handles.rawsignal - handles.bkg;
    % Replot
    switch handles.rawcorr
        case 'RAW'
            plot2D(handles,handles.rawsignal,handles.axes1,'On')
        case 'CORRECTED'
            plot2D(handles,handles.corrdata,handles.axes1,'On')
    end
    guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function maxtimeBkg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxtimeBkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotKin.
function plotKin_Callback(hObject, eventdata, handles)
interactive = get(handles.InteractiveModeTick,'Value');
handles.SelTraces = [];
switch interactive
    case 1
        % Select traces interactively, to finish hit RETURN
        handles = SelectTraces(handles);     
        % Check if RAW or CORR and give DATA
        switch handles.rawcorr
            case 'RAW'
                data = handles.rawsignal;
            case 'CORRECTED'
                data = handles.corrdata;
        end
    case 0
        handles.SelTraces = inputdlg('Enter space-separated wavenumbers to plot:',...
             'Input desired wavenumbers', [1 50]);
        handles.SelTraces = transpose(str2num(handles.SelTraces{:}));
        handles.SelTraces = sort(handles.SelTraces,1); 
        % Check if RAW or CORR and give DATA
        switch handles.rawcorr
            case 'RAW'
                data = handles.rawsignal;
            case 'CORRECTED'
                data = handles.corrdata;
        end
end
% ONCE EVERYTHING IS THERE, THEN DO STUFF
% Get the values to find in the WAVELENGTH vector
value = handles.SelTraces(:,1); % column 1= Wavelength
% k = index of the closest match for each selected trace
k = unique(findClosestId2Val(handles.cmprobe,transpose(value)));
L = length(k);
i = 0;
while i < L
    % Get the actual vectors from the Z matrix (one by one)
    index = k(1,i+1);
    switch handles.NormKin.Value
        case 0
            ydata(:,i+1) = data(:,index);
            label = '\DeltaAbs (mOD)';
            caption(i+1,1) = string(strcat(num2str(round(handles.cmprobe(index),0)),' cm^{-1}'));
        case 1
            maxval = max(data(:,index));
            minval = min(data(:,index));
            if abs(maxval) >= abs(minval)
                ydata(:,i+1) = data(:,index)./maxval;
                caption(i+1,1) = string(strcat(num2str(round(handles.cmprobe(index),0)),' cm^{-1}'));
            else
                ydata(:,i+1) = data(:,index)./minval;
                caption(i+1,1) = string(strcat(num2str(round(handles.cmprobe(index),0)),' cm^{-1} \times -1'));
            end
            label = 'Normalised \DeltaAbs';
    end
    i=i+1;
end
% Create a new figure with consistent format
fh = figure();
fh.Position(3)  = 920;
fh.Position(4)  = 425;
fh.Color        = [1 1 1];

% Define the axes
handles.axes2 = axes('Parent',fh);
axes(handles.axes2);

% % Plot in two axes 
% figure;
% ax1 = axes('OuterPosition',[0 0 0.3 1],'Units','Normalized');
% ax2 = axes('OuterPosition',[0.1 0 0.8 1],'Units','Normalized','YTickLabel',[]);

% Plot the data
cmap=colormap((othercolor('Mrainbow',L)));
for n=1:L
   plot(handles.delays,ydata(:,n),'LineWidth',2,'Marker','o','MarkerSize',2,'color',cmap(n,:));
   hold on
end
% Nice formatting
set(gca,'FontSize',12)
title({handles.datafilename;[handles.rawcorr,' DATA - KINETICS';'']},'Interpreter','none','FontSize',12)
xlabel(['Delays',' (',handles.timescale,')'],'FontSize',13,'FontWeight','bold');
ylabel(label,'FontSize',13,'FontWeight','bold')
axis tight
if L == 1
    legend(gca,[string(strcat(num2str(round(handles.cmprobe(index),0)),' cm^{-1}')),''],'FontSize',12)
else
    legend(gca,caption,'FontSize',12)
end
legend('boxoff')
legend('Location','bestoutside')
set(handles.axes2,'xscale',handles.linlog)
hline = refline(handles.axes2,[0 0]); hline.Color = [0.5 0.5 0.5];
set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
guidata(hObject,handles)
handles.axes2.Units     = 'pixels';
handles.axes2.Position  = [75 60 675 320];
handles.axes2.Units     = 'normalized';

% Decide whether to save the plotted traces or not
if handles.DoSaveTraces==1
    wavenumbers=transpose(handles.cmprobe(k));
    filename=char(strcat(handles.CurrDir.String,filesep,handles.datafilename,'_traces.dat'));
    dlmwrite(filename,[[0;handles.delays] [wavenumbers;ydata]])
end

% --- Executes on button press in plotDeltaAbs.
function plotDeltaAbs_Callback(hObject, eventdata, handles)
% hObject    handle to plotKin (see GCBO)-
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
interactive = get(handles.InteractiveModeTick,'Value');
handles.SelTraces = [];
switch interactive
    case 1
        timescale = handles.timescale;
        % Select traces interactively, to finish hit RETURN
        handles = SelectTraces(handles);     
        % Check if RAW or CORR and give DATA
        switch handles.rawcorr
            case 'RAW'
                data = transpose(handles.rawsignal);
            case 'CORRECTED'
                data = transpose(handles.corrdata);
        end
        % Get the values to find in the TIME vector
        value = handles.SelTraces(:,2); % column 2 = Time
    case 0
        timescale = handles.timescale;
        handles.SelTraces = inputdlg('Enter space-separated times to plot the transient spectra:',...
             ['Input desired times (in ',timescale,' )'], [1 50]);
        handles.SelTraces = transpose(str2num(handles.SelTraces{:}));
        handles.SelTraces = sort(handles.SelTraces,1); 
        % Check if RAW or CORR and give DATA
        switch handles.rawcorr
            case 'RAW'
                data = transpose(handles.rawsignal);
            case 'CORRECTED'
                data = transpose(handles.corrdata);
        end
         % Get the values to find in the TIME vector
        value = handles.SelTraces(:,1); % column 1 = Time (the only one)
end
% ONCE EVERYTHING IS THERE, THEN DO STUFF
% k = index of the closest match for each selected trace
k = findClosestId2Val(handles.delays,transpose(value));
k = unique(k,'stable');
L = length(k);
% ydata   = zeros(size
for n=1:L
    % Check the delays and convert the labels and numbers to appropriate units
    if handles.delays(k(n)) >= 1000
       newdelays = handles.delays(k(n))/1000;
       switch handles.timescale
           case 'fs'
               newtimescale = 'ps';
           case 'ps'
               newtimescale = 'ns';
           case 'ns'
               newtimescale = ['\mu' 's'];
       end
    elseif handles.delays(k(n)) < 1
       newdelays = handles.delays(k(n)).*1000;
       switch handles.timescale
           case 'ps'
               newtimescale = 'fs';
           case 'ns'
               newtimescale = 'ps';
       end
    else
       newdelays    = handles.delays(k(n));
       newtimescale = handles.timescale;
    end
    
    % Get the actual vectors from the Z matrix (one by one)
    % together with their legend captions and corresponding ylabel
    switch handles.NormKin.Value
        case 0
            ydata(:,n) = data(:,k(n));
            label = '\DeltaAbs (mOD)';
        case 1
            ydata(:,n) = data(:,k(n))./max(max(abs(data(:,k(n)))));
            label = 'Normalised \DeltaAbs';
    end
    caption(n,1) = string([num2str(round(newdelays,2,'significant')) ' ' newtimescale]);
end
   
switch handles.IncludeSteadyState.Value
    case 1
        % Create a new figure with consistent format
        fh = figure();
        fh.Units        = 'pixels';
        fh.Position(2)  = fh.Position(2)-150;
        fh.Position(3)  = 880;
        fh.Position(4)  = 550;
        fh.Color        = [1 1 1];
        % Create the FTIR axis
        ax1             = axes('Parent',fh);
        ax1.Units       = 'pixels';
        ax1.Position    = [75 400 675 125];
        plot(ax1,handles.IRx,handles.IRy,'LineWidth',1.5,'Marker','none','color','black');
%         title(ax1,{handles.datafilename;[handles.rawcorr,' DATA - TRANSIENT SPECTRA';'']},'Interpreter','none')
        ylabel(ax1,'Abs (mOD)','FontWeight','bold','FontSize',13)
        ax1.Units       = 'normalized';
        ax1.FontSize    = 12;
        % Create the transient spectra axis
        ax2             = axes('Parent',fh);
        ax2.FontSize    = 12;
        linkaxes([ax1,ax2],'x');
        where = ax2;
    case 0
        % Create a new figure with consistent format
        fh = figure();
        fh.Position(3)  = 880;
        fh.Position(4)  = 425;
        fh.Color        = [1 1 1];
        % Define the axes
        ax2 = axes('Parent',fh);
        axes(ax2);
        ax2.Units     = 'pixels';
        title(ax2,{handles.datafilename;[handles.rawcorr,' DATA - TRANSIENT SPECTRA';'']},'Interpreter','none')
        where = ax2;
end
% Plot the data
cm=colormap(othercolor('Mrainbow',L));
for n=1:L
   plot(where,handles.cmprobe,ydata(:,n),'LineWidth',1.5,'Marker','none','color',cm(n,:),'DisplayName',caption(n));
   hold on
end
% Nice formatting
set(gca,'FontSize',12)
xlabel('Wavenumbers (cm^{-1})','FontSize',13,'FontWeight','bold')
ylabel(label,'FontSize',13,'FontWeight','bold')
axis tight
if handles.IncludeSteadyState.Value == 1
    xl = xlim;
    xlim(ax1,xl);
    set(ax1, 'XTick', []);
end
legend('show');
legend('boxoff')
legend('Location','bestoutside')

where.Units       = 'pixels';
where.Position    = [75 75 675 320];
where.Units       = 'normalized';

if L >= 15
    for i=1:floor(L/3)
        where.Children(2*i).Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end

legend(ax2,{},'FontSize',12)

% Disable warning of negative X limits
warning('off','MATLAB:Axes:NegativeLimitsInLogAxis');

% Refline
hline = refline(0,0); hline.Color = [0.5 0.5 0.5];
set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
guidata(hObject,handles)

% --- Executes on button press in ExitButton.
function ExitButton_Callback(hObject, eventdata, handles)
% hObject    handle to ExitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all


% --- Executes on button press in ShowContoursTick.
function ShowContoursTick_Callback(hObject, eventdata, handles)
switch handles.rawcorr
    case 'RAW'
        plot2D(handles,handles.rawsignal,handles.axes1,'On')
    case 'CORRECTED'
        plot2D(handles,handles.corrdata,handles.axes1,'On')
end
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of ShowContoursTick


% --- Executes on button press in AddReplace1Dplot.
function AddReplace1Dplot_Callback(hObject, eventdata, handles)
% hObject    handle to AddReplace1Dplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(hObject,'Value')
    case 0
        % REPLACE
    case 1
        % ADD
end
guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of AddReplace1Dplot


% --- Executes on button press in BkgSubTick.
function BkgSubTick_Callback(hObject, eventdata, handles)
% hObject    handle to BkgSubTick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
switch value
    case 1
        handles.rawcorr = 'CORRECTED';
        plot2D(handles,handles.corrdata,handles.axes1,'On')
    case 0
        handles.rawcorr = 'RAW';
        plot2D(handles,handles.rawsignal,handles.axes1,'On')
end
guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of BkgSubTick



function ContourLineStyle_Callback(hObject, eventdata, handles)
switch handles.rawcorr
    case 'RAW'
        plot2D(handles,handles.rawsignal,handles.axes1,'On')
    case 'CORRECTED'
        plot2D(handles,handles.corrdata,handles.axes1,'On')
end
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ContourLineStyle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ContourLineStyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlotNoise.
function PlotNoise_Callback(hObject, eventdata, handles)
% hObject    handle to PlotNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
linlog = handles.linlog;
fh              = figure();
newaxes         = axes('Parent',fh);
plot3D(handles,handles.noise,newaxes)
set(newaxes,'yscale',linlog)

% Make the plot format constant
fh.Units='pixels';
fh.Position(2)=fh.Position(2).*0.7;
fh.Position(3)=620; % Width in pixels
fh.Position(4)=510; % Height in pixels
fh.Color = [1 1 1];
newaxes.FontSize = 12;
newaxes.Units = 'pixels';
newaxes.Position = [75 75 420 385];
% Make uniform, consistent format
newaxes.LineWidth = 1;
newaxes.TickLength = [0.015 0.035];
box(newaxes,'on')
% Make the axes resizable
newaxes.Units = 'normalized';

% Customise the graph a bit
maxnoise = abs(max(handles.noise));
minnoise =  abs(min(handles.noise));
Zminmax = max([maxnoise minnoise]);
caxis([0,Zminmax])
title({handles.datafilename;' NOISE'},'Interpreter','none')
view(0,-90)
colormap(othercolor('Spectral10'));
guidata(hObject,handles)


% --- Executes on selection change in DataTypeMenu.
function DataTypeMenu_Callback(hObject, eventdata, handles)
% hObject    handle to DataTypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns DataTypeMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DataTypeMenu


% --- Executes during object creation, after setting all properties.
function DataTypeMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataTypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in KineticsPerScan.
function KineticsPerScan_Callback(hObject, eventdata, handles)

Nscans          = handles.Nscans;
datafilename    = char(handles.datafilename);
tempdir         = [char(handles.CurrDir.String) filesep datafilename filesep 'temp'];
tempdirOUT      = [char(handles.CurrDir.String) filesep datafilename 'temp'];
% Check if the temp dir is outside (Lab 1)
if exist(tempdir,'dir') == 0
    tempdir     = tempdirOUT;
end

cd (char(tempdir));
interactive = get(handles.InteractiveModeTick,'Value');
handles.SelTraces = [];
switch interactive
    case 1
        % Select ONE point to compare the KINETICS
        hold on
        [xp,yp,button] = ginput(1);
        if isequal(button,1)
            plot(xp,yp,'X','Linewidth',1.5,'Color','k')
            handles.SelTraces = xp;
        end
        hold off
    case 0
        handles.SelTraces = inputdlg('Enter ONE wavenumber to plot:',...
             'Input desired wavenumber', [1 20]);
        handles.SelTraces = str2num(char(handles.SelTraces));
end
% Once we have the point, look for the index
% Get the values to find in the WAVELENGTH vector
value = handles.SelTraces; % column 1= Wavelength = chosen x coordinate
% k = index of the closest match for the selected trace
k = findClosestId2Val(handles.cmprobe,value);
% Load the scan data depending on Lab 1 vs Lab 2
switch handles.DataTypeMenu.Value
    case 1 % Lab 2
        % The first scan has no _n.csv termination
            Scan1 = char(strcat(tempdir,filesep,datafilename,'_signal.csv'));
            TempScanData = csvread(Scan1);
        if Nscans > 1
            ScanData(:,1) = TempScanData(:,k);
            % Load scans 2 to Nscans
            for i=2:Nscans
                ScanNr = strcat(tempdir,filesep,datafilename,'_signal_',num2str(i-1),'.csv');
                ScanNrch = char(ScanNr);
                TempScanData = csvread(ScanNrch);
                ScanData(:,i) = TempScanData(:,k);
                caption(i) = string(['Scan ' num2str(i)]);   
            end
            caption(1)= "Scan 1";
            clear TempScanData;
        else
            caption = "Scan 1";
            ScanData(:,1) = TempScanData(:,k);
        end
    case 2 % Lab 1 & Lab 4
        % The first scan has no _n.csv termination
        ending = '_sp0_sm0_du0_';
        for i=1:Nscans
            ScanNr   = strcat(tempdir,filesep,datafilename,'_signal',ending,num2str(i-1),'.csv');
            ScanNrch = char(ScanNr);
            TempScanData = csvread(ScanNrch);
            % Average rows 1:n and subtract them
            n=3; 
            TempScanData =  TempScanData(n:end,:) - mean(TempScanData(1:n-1,:),1);
            if handles.removeduplicates == 1
%                 % Average last 3 rows and collapse into a single row
%                 TempScanData(end-2,:) = mean(TempScanData(end-2:end,:),1);
%                 TempScanData(end-1:end,:)=[];
            end
            % Store it
            ScanData(:,i) = TempScanData(:,k);
            % Caption it
            caption(i) = string(['Scan ' num2str(i)]);   
        end
        clear TempScanData;
end
Nplots = Nscans;
% Bin the scans if BinScans > 1
BinScans = str2double(handles.BinScans.String);
if BinScans > 1
    s=1; i=1;
    BinScanData=zeros(length(handles.delays),ceil(Nscans/BinScans));
    while s < Nscans
        snew = s+BinScans;
        if snew > Nscans 
            m=[s:Nscans];
            BinCaption(i)=string([num2str(s) ' to ' num2str(Nscans)]);
        else
            m=[s:snew-1];
            BinCaption(i)=string([num2str(s) ' to ' num2str(snew-1)]);
        end
        BinScanData(:,i)=mean(ScanData(:,m),2);
        s=snew; i=i+1;
    end
    RawScanData = ScanData;
    ScanData    = BinScanData;
    Nplots      = ceil(Nscans/BinScans);
    caption     = BinCaption;
end

% Normalise the data if NormKin = 1
switch handles.NormKin.Value
    case 0
        label = '\DeltaAbs (mOD)';
    case 1
        label = 'Normalised \DeltaAbs';
        j=0;
        while j < Nplots
            j=j+1;
            ScanData(:,j) = ScanData(:,j)./max(max(abs(ScanData(:,j))));
        end
end

% Create a new figure
fh = figure();
fh.Position(3)  = 920;
fh.Position(4)  = 425;
fh.Color        = [1 1 1];

% Define the axes
handles.axes2 = axes('Parent',fh);
axes(handles.axes2);

% Plot the data
cm=colormap(othercolor('Mrainbow',Nplots));
for n=1:Nplots
   plot(handles.delays,ScanData(:,n),'color',cm(n,:));
   hold on
end
% Make the lines of the first and last scans thicker
handles.axes2.Children(end).LineWidth = 1.5;
handles.axes2.Children(1).LineWidth = 1.5;

% Nice formatting
title({handles.datafilename;['KINETICS PER SCAN AT ',num2str(round(handles.cmprobe(k),0)),' cm-1']},'Interpreter','none')
xlabel(['Delays',' (',handles.timescale,')'],'FontSize',13,'FontWeight','bold')
ylabel(label,'FontSize',13,'FontWeight','bold')
axis tight
legend(gca,caption);
legend(gca,'Location','bestoutside');
legend(gca,'boxoff');
set(handles.axes2,'xscale',handles.linlog)
hline = refline(handles.axes2,[0 0]);
hline.Color = [0.5 0.5 0.5];
set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
handles.axes2.FontSize  = 12;
handles.axes2.Units     = 'pixels';
handles.axes2.Position  = [75 60 675 320];
handles.axes2.Units     = 'normalized';
guidata(hObject,handles)


% --- Executes on button press in SpectraPerScan.
function SpectraPerScan_Callback(hObject, eventdata, handles)
% hObject    handle to KineticsPerScan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datafilename    = char(handles.datafilename);
tempdir         = [char(handles.CurrDir.String) filesep datafilename filesep 'temp'];
tempdirOUT      = [char(handles.CurrDir.String) filesep datafilename 'temp'];
% Check if the temp dir is outside (Lab 1)
if exist(tempdir,'dir') == 0
    tempdir     = tempdirOUT;
end
cd (char(tempdir));
filelist = dir('*signal*.csv');
Nscans = length(filelist)./2;
interactive = get(handles.InteractiveModeTick,'Value');
handles.SelTraces = [];
switch interactive
    case 1
        % Select ONE point to compare the KINETICS
        hold on
        [xp,yp,button] = ginput(1);
        if isequal(button,1)
            plot(xp,yp,'X','Linewidth',1.5,'Color','k')
            handles.SelTraces = yp;
        end
        hold off
    case 0
        handles.SelTraces = inputdlg('Enter ONE delay to plot:',...
             'Input desired wavenumber', [1 20]);
        handles.SelTraces = str2num(char(handles.SelTraces));
end
% Once we have the point, look for the index
% Get the values to find in the TIME vector
value = handles.SelTraces; % column 1= Wavelength = chosen x coordinate
% k = index of the closest match for the selected trace
k = findClosestId2Val(handles.delays,value);
switch handles.DataTypeMenu.Value
    case 1 % Lab 2
        % The first scan has no _n.csv termination
            Scan1 = char(strcat(tempdir,filesep,datafilename,'_signal.csv'));
            TempScanData = csvread(Scan1);
        if Nscans > 1
            ScanData(:,1) = TempScanData(k,:);
            % Load scans 2 to Nscans
            for i=2:Nscans
                ScanNr = strcat(tempdir,filesep,datafilename,'_signal_',num2str(i-1),'.csv');
                ScanNrch = char(ScanNr);
                TempScanData = csvread(ScanNrch);
                ScanData(:,i) = TempScanData(k,:);
                caption(i) = string(['Scan ' num2str(i)]);   
            end
            caption(1)= "Scan 1";
            clear TempScanData;
        else
            caption = "Scan 1";
            ScanData(:,1) = TempScanData(:,k);
        end
        ScanData = transpose(ScanData);
    case 2 % Lab 1 & Lab 4
        % The first scan has no _n.csv termination
        ending = '_sp0_sm0_du0_';
        Nscans = handles.Nscans;
        for i=1:Nscans
            ScanNr   = strcat(tempdir,filesep,datafilename,'_signal',ending,num2str(i-1),'.csv');
            ScanNrch = char(ScanNr);
            TempScanData = csvread(ScanNrch);
            if handles.removeduplicates == 1
                % Average rows 1:n and subtract them (n=5)
                n=5; 
                TempScanData =  TempScanData(n:end,:) - mean(TempScanData(1:n-1,:),1);
                % Average last 3 rows and collapse into a single row
                TempScanData(end-2,:) = mean(TempScanData(end-2:end,:),1);
                TempScanData(end-1:end,:)=[];
            end
            % Store it
            ScanData(i,:) = TempScanData(k,:);
            % Caption it
            caption(i) = string(['Scan ' num2str(i)]);   
        end
        clear TempScanData;
end
Nplots = Nscans;

% Bin the scans if BinScans > 1
BinScans = str2double(handles.BinScans.String);
if BinScans > 1
    s=1; i=1;
    BinScanData=zeros(ceil(Nscans/BinScans),length(handles.cmprobe));
    while s <= Nscans
        snew = s+BinScans;
        if snew > Nscans 
            m=s:Nscans;
            BinCaption(i)=string([num2str(s) ' to ' num2str(Nscans)]);
        else
            m=s:snew-1;
            BinCaption(i)=string([num2str(s) ' to ' num2str(snew-1)]);
        end
        BinScanData(i,:)=mean(ScanData(m,:),1);
        s=snew; i=i+1;
    end
    RawScanData = ScanData;
    ScanData = BinScanData;
    Nplots = ceil(Nscans/BinScans);
    caption = BinCaption;
end

% Normalise the data if NormKin = 1
switch handles.NormKin.Value
    case 0
        label = '\DeltaAbs (mOD)';
    case 1
        label = 'Normalised \DeltaAbs';
        for j=1:Nplots
            ScanData(j,:) = ScanData(j,:)./max(max(abs(ScanData(j,:))));
        end
end

% Create a new figure with consistent format
fh = figure();
fh.Position(3)  = 880;
fh.Position(4)  = 425;
fh.Color        = [1 1 1];

% Define the axes
handles.axes2 = axes('Parent',fh);
axes(handles.axes2);

% Plot the data
cm=colormap(othercolor('Mrainbow',Nplots));
for n=1:Nplots
   pl(n) = plot(handles.cmprobe,ScanData(n,:),'color',cm(n,:));
   hold on
end
% Make the lines of the first and last scans thicker
pl(end).LineWidth   = 1.5;
pl(1).LineWidth     = 1.5;

% Nice formatting
set(handles.axes2,'FontSize',12)
title({handles.datafilename;['TR. SPECTRA PER SCAN AT ',num2str(handles.delays(k),'%.3g'),' ' handles.timescale]},'Interpreter','none','FontSize',12)
xlabel('Wavenumbers (cm^{-1})','FontSize',13,'FontWeight','bold')
ylabel(label,'FontSize',13,'FontWeight','bold')
axis tight
hline       = refline(0,0); 
hline.Color = [0.5 0.5 0.5];
legend(gca,caption);
legend(gca,'boxoff');
legend(gca,'Location','bestoutside');
handles.axes2.Units     = 'pixels';
handles.axes2.Position  = [75 60 675 320];
handles.axes2.Units     = 'normalized';
guidata(hObject,handles)


% --- Executes on button press in InteractiveModeTick.
function InteractiveModeTick_Callback(hObject, eventdata, handles)
% THIS TICK SHOULD NOT DO ANYTHING

% --- Executes on button press in NormKin.
function NormKin_Callback(hObject, eventdata, handles)
% THIS TICK SHOULD NOT DO ANYTHING

% --- Executes on button press in Refresh.
function Refresh_Callback(hObject, eventdata, handles)
rootdir = handles.CurrDir.String;

switch handles.DataTypeMenu.Value
    case 4 % *.dat files
        if rootdir ~= 0
            % Update text accordingly
            handles.CurrDir.String  = rootdir;
            % Get a list of all files and folders in this folder.
            files = dir(rootdir);
            % Get only those that are *.dat
            files = files(contains({files.name},".dat",'IgnoreCase',1));
            % Sort the names and create variables for the list
                [sorted_names,sorted_index] = sortrows({files.name}');
                handles.file_names          = sorted_names;
                handles.sorted_index        = sorted_index;
                guidata(hObject,handles)
                set(handles.DatafoldersList,'String',handles.file_names,'Value',1)
                handles.is_dir = 0;
                guidata(hObject,handles)
        else 
            handles.CurrDir.String = "Please select a directory!  :(";
        end
    otherwise % = 1,2,3
        if rootdir ~= 0
            % Update text accordingly
            handles.CurrDir.String  = rootdir;
            % Get a list of all files and folders in this folder.
            files = dir(rootdir);
            % Extract only those that are directories.
            subFolders = files([files.isdir]);
%           % Go to the folder
%           cd (rootdir);
            % Sort the names and create variables for the list
                [sorted_names,sorted_index] = sortrows({subFolders.name}');
                handles.file_names = sorted_names;
                handles.is_dir = [subFolders.isdir];
                handles.sorted_index = sorted_index;
                guidata(hObject,handles)
            % Update DatafoldersList, removing ".." and "." from the visible list
                handles.file_names(1:2,:) = [];
                set(handles.DatafoldersList,'String',handles.file_names,'Value',1)
                guidata(hObject,handles)     
        %       else
        %         % LOAD DATA IN LAB 4 PUMP-PROBE FORMAT
        %         % Execute PatchFileList_Lab4 and use the output to build the file list
        %             handles = PatchFileList_Lab4(handles,'1D');
        %             guidata(hObject,handles);
        %             set(handles.DatafoldersList,'String',handles.SampleName,'Value',1);
        %       end
        else 
            set(handles.CurrDir, 'String', 'Please select a directory!  :(');
        end
end

% --- Executes on button press in Plot_Dependency.
function Plot_Dependency_Callback(hObject, eventdata, handles)
% Ask the user what to plot
plot_options = {'Concentration dependence','Power dependence'};
[plot_typeindx,doplot] = listdlg('ListString',plot_options,'OKstring','Plot','SelectionMode','single','ListSize',[180,60],'PromptString','Select slice type to plot:');

if doplot == 0
    return
end

switch plot_options{plot_typeindx}
    case 'Concentration dependence'
        % Get the data from user - The user should put all data in one directory,
        % then indicate the concentration for each dataset
        AllFolderList = handles.DatafoldersList.String;
        rootdir = handles.CurrDir.String;
        % Ask user for the wavenumber to plot and the times to plot it
        prompt = {'Enter desired wavenumber(s) to plot, separated by a space:',['Enter time range: (in ' handles.timescale ')']};
        dlg_title = 'Concentration dependence of signal';
        num_lines = [1 30];
        defAns = {'',[num2str(handles.delays(1)) ' ' num2str(handles.delays(end))]};
        SelTraces = inputdlg(prompt,dlg_title,num_lines,defAns);
        WL = str2num(SelTraces{1});
        Times = str2num(SelTraces{2});
        % Once we have the point, look for the index
        % WLindex = index of the closest match for the selected trace in the WL vector
        WLindex = findClosestId2Val(handles.cmprobe,WL);
        timeindex = findClosestId2Val(handles.delays,Times);
        tmin = timeindex(1);
        tmax = timeindex(2);
        % Then, ask the user to input the concentration for each trace in the current rootdir.
        % To skip a particular trace, enter -1. The conc. should be typed in the
        % same order as the FoldersList is displayed
        prompt = {'Enter the concentration values in the same order as the file list (to skip a datafile use -1):','Enter concentration units:'};
        dlg_title = 'Plot of conc. dependence of signal';
        num_lines = [1 50];
        DlgConc = inputdlg(prompt,dlg_title,num_lines,{'','mM'});
        Concentrations = transpose(str2num(DlgConc{1}));
        ConcUnits = DlgConc{2};
        Nconc = length(Concentrations);
        SelConc=[];
        % Remove all folders which have a -1 from the list
        i=0; j=1;
        while i < Nconc
            i=i+1;
            if Concentrations(i) ~= -1
                SelConc(j) = Concentrations(i);
                FolderList(j,1) = AllFolderList(i);
                j=j+1;
            end
        end    
        Nconc = length(SelConc);
        [~,I]=sort(SelConc);
        % Get the actual vectors from the Z matrix (one by one) and plot them all
        % together for a given WL
        for WLn=1:length(WL)
            fh = figure('Units','normalized','Position',[0.25+WLn*0.05 0.25+WLn*0.05 0.35 0.35]);
            if WLn==1
                cm=colormap(othercolor('Mrainbow',Nconc));
            end
            for p=1:Nconc
                    datafilename = FolderList(I(p));
                switch handles.DataTypeMenu.Value
                    case 1 % Lab 2
                        signalfile = strcat(rootdir,filesep,datafilename,filesep,datafilename,'_signal.csv');
                    case 2 % Lab 1 and Lab 4
                        signalfile = strcat(rootdir,filesep,datafilename,filesep,datafilename,'_signal_sp0_sm0_du0.csv');
                end
                    signalfilech = char(signalfile);
                    tempsignal = csvread(signalfilech);
                    % Auto subtract background
                    k = 1;
                    t = handles.delays(k);
                    while t < -1
                        k = k+1;
                        t = handles.delays(k+1);
                    end
                    handles.j = 1;
                    handles.k = k;
                    set(handles.mintimeBkg,'String',num2str(handles.delays(1)));
                    set(handles.maxtimeBkg,'String',num2str(handles.delays(k)));
                    % Do the background subtraction
                    bkg = mean(tempsignal(1:k,:));
                    corrdata = tempsignal - bkg;
                    trace=corrdata(tmin:tmax,WLindex(WLn));
                    if handles.NormKin.Value == 1
                        trace=trace/max([abs(min(trace(:))),abs(max(trace(:)))]);
                    end
                    % Plot the traces as a function of concentration for the current WL
                    hold on
                    plot(handles.delays(tmin:tmax),trace,'Color',cm(p,:),'LineWidth',1.5,'DisplayName',[num2str(SelConc(I(p))) ' ' ConcUnits])
            end
            nicetitle = {['CONCENTRATION DEPENDENCE PLOT - SIGNAL AT ',num2str(handles.cmprobe(WLindex(WLn)),'%4.f'),' cm^-^1'];['from ',num2str(handles.delays(tmin)),' to ',num2str(handles.delays(tmax)),' ',handles.timescale]}; 
            title(nicetitle,'FontSize',10)
            set(gca,'fontsize',12);
            xlabel(['Delays (' handles.timescale ')'],'FontWeight','bold','FontSize',13)
            if handles.NormKin.Value == 1
                ylabel('Normalized \DeltaAbs','FontWeight','bold','FontSize',13)
            else
                ylabel('\DeltaAbs (mOD)','FontWeight','bold','FontSize',13)
            end
            axis tight
            legend(gca,'Location','best');
            legend('boxoff');
            set(gca,'xscale',handles.linlog);
            box(gca,'on');
            limits=xlim(gca);
            % Add reference lines at zero
            hline=refline(0,0);
            hline.Color = [0.5 0.5 0.5];
            hold(gca,'off')
            ax = gca;
            ax.Legend.String(end)=[];

            % Make the plot format constant
            % Create a new figure with consistent format
            fh.Units        = 'pixels';
            fh.Position(3)  = 800;
            fh.Position(4)  = 425;
            fh.Color        = [1 1 1];

            % Make uniform, consistent format
            ax.Units     = 'pixels';
            ax.Position  = [75 55 680 320];  
            ax.FontSize = 12;
            ax.LineWidth = 1;
            ax.TickLength = [0.015 0.035];

            % Make the axes resizable
            ax.Units = 'normalized';   
        end
    case 'Power dependence'
        % Get the data from user - The user should put all data in one directory, then indicate the power for each measurement
        AllFolderList = handles.DatafoldersList.String;
        rootdir = handles.CurrDir.String;
        % Ask user for the wavenumber to plot and the times to plot it
        prompt      = {'Enter desired wavenumber:','Enter time range to average (min):','Enter time range to average (max):'};
        dlg_title   = 'Plot of power dependence of signal';
        num_lines   = [1 20];
        SelTraces   = inputdlg(prompt,dlg_title,num_lines);
        WL          = str2double(SelTraces(1));
        Times       = transpose(str2double(SelTraces(2:3)));
        % Once we have the point, look for the index
        % WLindex = index of the closest match for the selected trace in the WL vector
        WLindex = findClosestId2Val(handles.cmprobe,WL);
        timeindex = findClosestId2Val(handles.delays,Times);
        tmin = timeindex(1);
        tmax = timeindex(2);
        % Then, ask the user to input the power for each trace in the current rootdir.
        % To skip a particular trace, enter -1. The power should be typed in the
        % same order as the FoldersList is displayed
        prompt = {'Enter the power values in the same order as the list (to skip a datafile use -1):'};
        dlg_title = 'Plot of power dependence of signal';
        num_lines = [1 50];
        DlgPowers = inputdlg(prompt,dlg_title,num_lines);
        DlgPowers = transpose(str2num(DlgPowers{:}));
        Npower = length(DlgPowers);
        % Remove all folders which have a -1 from the list
        i=0; j=1;
        while i < Npower
            i=i+1;
            if DlgPowers(i) ~= -1
                SelPowers(j) = DlgPowers(i);
                FolderList(j,1) = AllFolderList(i);
                j=j+1;
            end
        end    
        Npower = length(SelPowers);
        Avg=[];
        % Get the actual vectors from the Z matrix (one by one)
        p=0;
        while p < Npower
                p=p+1;
                datafilename = FolderList(p);
                signalfile = strcat(rootdir,filesep,datafilename,filesep,datafilename,'_signal.csv');
                signalfilech = char(signalfile);
                tempsignal = csvread(signalfilech);
                % Auto subtract background
                k = 1;
                t = handles.delays(k);
                while t < 0
                    k = k+1;
                    t = handles.delays(k+1);
                end
                handles.j = 1;
                handles.k = k;
                set(handles.mintimeBkg,'String',num2str(handles.delays(1)));
                set(handles.maxtimeBkg,'String',num2str(handles.delays(k)));
                % Do the background subtraction
                bkg = mean(tempsignal(1:k,:));
                corrdata = tempsignal - bkg;
                % Average the data within given times
                Avg(:,p) = mean(corrdata(tmin:tmax,WLindex));
        end
        fh = figure();
        Avg = transpose(Avg);
        plot(SelPowers,Avg,'Marker','o')
        nicetitle = {['POWER DEPENDENCE PLOT - SIGNAL AT ',num2str(handles.cmprobe(WLindex)),' cm^-^1'];['from ',num2str(tmin),' to ',num2str(tmax),' ',handles.timescale]}; 
        title(nicetitle,'FontSize',10)
        set(gca,'fontsize',12);
        xlabel(['Energy (' '\mu' 'J)'],'FontSize',13)
        ylabel('\DeltaAbs (mOD)','FontSize',13)
        axis tight
end


% --- Executes on button press in Fit.
function Fit_Callback(hObject, eventdata, handles)
% Ask the user for the wavenumber(s) to fit and the time interval to fit
switch handles.InteractiveModeTick.Value
    case 0
        prompt = {'Enter desired wavenumber(s), separated by a space:','Enter min and max time to fit (default All):'};
        dlg_title = 'Fitting range';
        num_lines = [1 55];
        SelTraces = inputdlg(prompt,dlg_title,num_lines);
        WL = transpose(str2num(SelTraces{1}));
        Times = transpose(str2num(SelTraces{2}));
    case 1
        % Select traces interactively, to finish hit RETURN
        handles = SelectTraces(handles);
        SelTraces = handles.SelTraces;
        WL = SelTraces(:,1);
        Times=[];
end
% Once we have the point, look for the index
% WLindex = index of the closest match for the selected trace in the WL vector
% timeindex = index of the closest match for the selected time in the delays vector
WLindex = findClosestId2Val(handles.cmprobe,WL);
% Tmin and Tmax are the indices of the tmin and tmax in the delays vector
if isempty(Times) % If the user gives no response in time, take all delays
    tmin = 1;
    tmax = length(handles.delays);
else
    timeindex = findClosestId2Val(handles.delays,Times);
    tmin = timeindex(1);
    tmax = timeindex(2);
end
% Determine whether the data is RAW or CORRECTED
switch handles.rawcorr
    case 'RAW'
        data=handles.rawsignal;
    case 'CORRECTED'
        data=handles.corrdata;
end
% Then, open the fit window
TraceWLs = handles.cmprobe(WLindex);
TraceLabels = num2str(TraceWLs,6);
data = data(:,WLindex);
YLabel = '\DeltaAbs (mOD)';
FitType=1; % i.e. 1= LOCAL, 2 = GLOBAL SVD fit
Output = FitWindow(handles.delays(tmin:tmax),data,TraceLabels,YLabel,handles.timescale,FitType);
% Save initial parameters to file
%formatSpec = 'Taus:\t\t[%6.2f\t%6.2f\t%6.2f\t%6.2f]\nTau_holds:\t\t[%i\t%i\t%i\t%i]\n Coeffs (1-4 inf):\t\t[%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t]\n Coeff_Holds(1-4 inf):\t\t[%i %i %i %i %i]\n FWHM:\t\t%6.2f\n t0:\t\t%6.2f\n Offset:\t\t%6.2f\n Holds(FWHM,t0,Offset):\t\t[%i %i %i]\n Do Fit?:\t\t%i\n';
%fprintf(formatSpec,handles.FitParam.Tau_values,handles.FitParam.Tau_holds,handles.FitParam.C_values,handles.FitParam.C_holds,handles.FitParam.Param_values,handles.FitParam.Param_holds,handles.FitParam.DoFit)


% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
m = length(handles.delays);
n = length(handles.cmprobe);
%The signal has m rows and n columns
XYZdata(1,2:n+1) = transpose(handles.cmprobe);
XYZdata(2:m+1,1) = handles.delays;
switch handles.rawcorr
    case 'RAW'
        XYZdata(2:m+1,2:n+1) = handles.rawsignal;
        datafilename = char(strcat(handles.datafilename,'-RAW'));
    case 'CORRECTED'
        XYZdata(2:m+1,2:n+1) = handles.corrdata;
        datafilename = char(strcat(handles.datafilename,'-CORRECTED'));
end
dlmwrite(strcat(datafilename,'.dat'),XYZdata);
helpdlg(['Data was saved in ' pwd filesep datafilename '.dat'],'Data saved successfully!');


% --- Executes on button press in UpDirectory.
function UpDirectory_Callback(hObject, eventdata, handles)
% hObject    handle to UpDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Go to the folder
cd (handles.CurrDir.String);
% Go one folder Up
cd ..;
% Update the string
handles.CurrDir.String = pwd;
rootdir = handles.CurrDir.String;
% Put list of subfolders into subFolders
    % Get a list of all files and folders in this folder.
    files = dir(rootdir);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
% Sort the names and create variables for the list
    dir_struct = subFolders;
    [sorted_names,sorted_index] = sortrows({dir_struct.name}');
    handles.file_names = sorted_names;
    handles.is_dir = [dir_struct.isdir];
    handles.sorted_index = [sorted_index];
    guidata(hObject,handles)
% Update DatafoldersList, removing ".." and "." from the visible list
    handles.file_names(1:2,:) = [];
    set(handles.DatafoldersList,'String',handles.file_names,'Value',1)
    guidata(hObject,handles)


% --- Executes on button press in LoadSpectrum.
function LoadSpectrum_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSpectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = LoadSteadyIR(handles);
guidata(hObject,handles)

% --- Executes on button press in tWLShift.
function tWLShift_Callback(hObject, eventdata, handles)
% hObject    handle to tWLShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Ask for the WL and t shifts
prompt = {'Enter the Wavelength and Time shifts (separated by a space, 0 to skip):'};
dlg_title = 'WL/t shift';
num_lines = [1 37];
DlgShift = inputdlg(prompt,dlg_title,num_lines);
DlgShift = transpose(str2num(DlgShift{:}));
WLShift = DlgShift(1);
tShift = DlgShift(2);
rootdir = handles.CurrDir.String;
% If not zero
if WLShift == 0 && tShift == 0
    set(handles.ErrorText,'String','Nothing Changed!');
else 
    % Write the values to handles
    handles.cmprobe = handles.cmprobe + WLShift;
    handles.delays = handles.delays + tShift;
    % Now, determine if the data comes from Lab4 (or not)
    if handles.DataTypeMenu.Value ~=3 % 1=lab2_pp, 2=lab1_pp, 3=lab4_pp, 4=lab2_truvis 
    % IF DATA DOESN'T COME FROM LAB 4
        % Update CSV's, creating a backup copy!
        datafilename = handles.datafilename;
        if WLShift ~= 0
            wavenumberfile = strcat(rootdir,filesep,datafilename,filesep,datafilename,'_wavenumbers.csv');
                newWNF = strcat(wavenumberfile,'bak');
                newWNFch = char(newWNF); wavenumberfilech = char(wavenumberfile);
                if exist(newWNFch,'file') == 0 %% Do not overwrite the .bak file!
                    movefile(wavenumberfilech,newWNFch);
                end
                csvwrite(wavenumberfilech,handles.cmprobe); %% Instead, overwrite only the wavenumberfile
                set(handles.ErrorText,'String','WL and t shifts written to disk!');
        end
        if tShift ~= 0
            delayfile = strcat(rootdir,filesep,datafilename,filesep,datafilename,'_delays.csv');
                newDF = strcat(delayfile,'bak');
                newDFch = char(newDF); delayfilech = char(delayfile);
                movefile(delayfilech,newDFch);
                csvwrite(delayfilech,handles.delays);
                set(handles.ErrorText,'String','WL and t shifts written to disk!');
        end
        handles = LoadDataIR2(handles);
    else
        % DATA COMES FROM LAB4
        ID = handles.DatafoldersList.Value;
        DelayZero = handles.DelayZero;
        if tShift ~= 0
            DelayZero(ID) = DelayZero(ID) - tShift;
            delayzerofile = [rootdir filesep 'DelayZero.txt'];
            csvwrite(delayzerofile,DelayZero)
            handles.DelayZero = DelayZero;
        end
    end
    switch handles.rawcorr
        case 'RAW'
            plot2D(handles,handles.rawsignal,handles.axes1,'On');
        case 'CORRECTED'
            plot2D(handles,handles.corrdata,handles.axes1,'On');
    end
end
% Update Handles
guidata(hObject,handles)


% --- Executes on button press in IncludeSteadyState.
function IncludeSteadyState_Callback(hObject, eventdata, handles)
% hObject    handle to IncludeSteadyState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
steady = handles.IncludeSteadyState.Value;
if steady == 1
    switch handles.FTIRloaded.Value
        case 0
            handles = LoadSteadyIR(handles);
            if handles.IRloaded == 0
              handles.IncludeSteadyState.Value = 0;
            end
    end
end
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of IncludeSteadyState


% --- Executes on button press in ShowFTIR.
function ShowFTIR_Callback(hObject, eventdata, handles)
% hObject    handle to ShowFTIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% New figure and axes
fh = figure();
handles.axes2 = axes('Parent',fh);
% Plot the data
plot(handles.IRx,handles.IRy,'LineWidth',1.5,'Marker','none','color','b');
% Nice formatting
set(gca,'FontSize',12)
title('FTIR SPECTRUM','FontSize',12)
xlabel('Wavenumbers (cm^{-1})','FontSize',13,'FontWeight','bold')
ylabel('Abs (mOD)','FontSize',13,'FontWeight','bold')
axis tight
hline = refline(0,0); hline.Color = [0.5 0.5 0.5];
guidata(hObject,handles)


% --- Executes on button press in EraseFTIR.
function EraseFTIR_Callback(hObject, eventdata, handles)
% hObject    handle to EraseFTIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear handles.IRx handles.IRy;
handles.FTIRloaded.Value = 0;
handles.IncludeSteadyState.Value = 0;
handles.ShowFTIR.Visible = 'Off';
handles.EraseFTIR.Visible = 'Off';
guidata(hObject,handles)

function nSVD_Callback(hObject, eventdata, handles)
% hObject    handle to nSVD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nSVD as text
%        str2double(get(hObject,'String')) returns contents of nSVD as a double

% --- Executes during object creation, after setting all properties.
function nSVD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nSVD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in SaveTraces.
function SaveTraces_Callback(hObject, eventdata, handles)
handles.DoSaveTraces = get(hObject,'Value');
guidata(hObject,handles)

% --- Executes on button press in PlotSumDiff.
function PlotSumDiff_Callback(hObject, eventdata, handles)
% hObject    handle to PlotSumDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
interactive = get(handles.InteractiveModeTick,'Value');
handles.SelTraces = [];
switch interactive
    case 1
        % Select traces interactively, to finish hit RETURN
        handles = SelectTraces(handles);     
        % Check if RAW or CORR and give DATA
        switch handles.rawcorr
            case 'RAW'
                data = handles.rawsignal;
            case 'CORRECTED'
                data = handles.corrdata;
        end
    case 0
        handles.SelTraces = inputdlg('Enter space-separated wavenumbers to plot:',...
             'Input desired wavenumbers', [1 50]);
        handles.SelTraces = transpose(str2num(handles.SelTraces{:}));
        % Check if RAW or CORR and give DATA
        switch handles.rawcorr
            case 'RAW'
                data = handles.rawsignal;
            case 'CORRECTED'
                data = handles.corrdata;
        end
end
% ONCE EVERYTHING IS THERE, THEN DO STUFF
% Get the values to find in the WAVELENGTH vector
value = handles.SelTraces(:,1); % column 1= Wavelength
% k = index of the closest match for each selected trace
k = findClosestId2Val(handles.cmprobe,transpose(value));
L = size(handles.SelTraces,1);
if L==2
    ydata = data(:,k);
    caption = string(strcat(num2str(round(handles.cmprobe(k),0)),' cm^{-1}'));
    
    % Create a new figure with consistent format
    fig1=figure();
    fig1.Position(3)  = 800;
    fig1.Position(4)  = 420;
    fig1.Color        = [1 1 1];
    % Define the axes
    newaxes1           = axes('Parent',fig1);
    axes(newaxes1);
    newaxes1.Units     = 'pixels';
    newaxes1.Position  = [75 60 680 320];
    box(newaxes1,'on');
    
    % % Plot in two axes 
    % figure;
    % ax1 = axes('OuterPosition',[0 0 0.3 1],'Units','Normalized');
    % ax2 = axes('OuterPosition',[0.1 0 0.8 1],'Units','Normalized','YTickLabel',[]);
    % Plot the data
    hold on
    % Difference: trace 1 - trace 2
    % Average of the whole spectrum
    Diff=ydata(:,1)-ydata(:,2);
    Avg=mean(data,2);
    Tr1_Avg = ydata(:,1)-Avg;
    Tr2_Avg = ydata(:,2)-Avg;
    % Normalise data if selected
    if handles.NormKin.Value == 1
        maxval_Diff = max(Diff);
        minval_Diff = min(Diff);
        maxval_Avg = max(Avg);
        minval_Avg = min(Avg);
        if abs(maxval_Diff) >= abs(minval_Diff)
            Diff = Diff./maxval_Diff;
        else
            Diff = Diff./minval_Diff;
        end
        if abs(maxval_Avg) >= abs(minval_Avg)
            Avg = Avg./maxval_Avg;
        else
            Avg = Avg./minval_Avg;
        end
        maxval_Tr1_Avg = max(Tr1_Avg);
        minval_Tr1_Avg = min(Tr1_Avg);
        maxval_Tr2_Avg = max(Tr2_Avg);
        minval_Tr2_Avg = min(Tr2_Avg);
        if abs(maxval_Tr1_Avg) >= abs(minval_Tr1_Avg)
            Tr1_Avg = Tr1_Avg./maxval_Tr1_Avg;
        else
            Tr1_Avg = Tr1_Avg./minval_Tr1_Avg;
        end
        if abs(maxval_Tr2_Avg) >= abs(minval_Tr2_Avg)
            Tr2_Avg = Tr2_Avg./maxval_Tr2_Avg;
        else
            Tr2_Avg = Tr2_Avg./minval_Tr2_Avg;
        end
        label = 'Normalised \DeltaAbs';
    else
        label = '\DeltaAbs (mOD)';
    end
    % Plot the data
    plot(handles.delays,Diff,'-or','LineWidth',2,'DisplayName',[caption{1,1} ' - ' caption{2,1}]);
    plot(handles.delays,Avg,'-xb','LineWidth',2,'MarkerSize',4,'DisplayName','Global average');
    plot(handles.delays,Tr1_Avg,'-*g','LineWidth',2,'MarkerSize',4,'DisplayName',[caption{1,1} ' - Global avg.']);
    plot(handles.delays,Tr2_Avg,'-*m','LineWidth',2,'MarkerSize',4,'DisplayName',[caption{2,1} ' - Global avg.']);
    hold off
    % Nice formatting
    set(gca,'FontSize',12)
    title({handles.datafilename;[handles.rawcorr,' DATA - KINETICS';'']},'Interpreter','none','FontSize',12)
    xlabel(['Delays',' (',handles.timescale,')'],'FontSize',13,'FontWeight','bold');
    ylabel(label,'FontSize',13,'FontWeight','bold')
    axis tight
    legend('show')
    legend('boxoff')
    legend('Location','best')
    set(newaxes1,'xscale',handles.linlog)
    title({handles.datafilename;[handles.rawcorr,' DATA';'']},'Interpreter','none')
    hline = refline(newaxes1,[0 0]); hline.Color = [0.5 0.5 0.5];
    set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    guidata(hObject,handles)
    
    % Plot the ratio between the molecular signals and the e- offset
    Ratio=Diff./Avg;
    % Create a new figure with consistent format
    fig2=figure();
    fig2.Position(3)  = 800;
    fig2.Position(4)  = 420;
    fig2.Color        = [1 1 1];
    % Define the axes
    newaxes2           = axes('Parent',fig2);
    axes(newaxes2);
    newaxes2.Units     = 'pixels';
    newaxes2.Position  = [75 60 680 320];
    box(newaxes2,'on');
    % Plot
    plot(newaxes2,handles.delays,Ratio,'-ok','LineWidth',2,'MarkerSize',4,'DisplayName','Ratio of Molec. signal to e- offset');
    % Prettify
    set(gca,'FontSize',12)
    title({handles.datafilename;[handles.rawcorr,' DATA - KINETICS';'']},'Interpreter','none','FontSize',12)
    xlabel(['Delays',' (',handles.timescale,')'],'FontSize',13,'FontWeight','bold');
    ylabel(label,'FontSize',13,'FontWeight','bold')
    axis tight
    set(gca,'XScale','log');
    legend('show')
    legend('boxoff')
    legend('Location','best')
    set(newaxes2,'xscale',handles.linlog)
    title({handles.datafilename;[handles.rawcorr,' DATA';'']},'Interpreter','none')
    hline = refline(newaxes2,[0 0]); hline.Color = [0.5 0.5 0.5];
    set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    guidata(hObject,handles)
else
    % Show a dialog that says that you must select 2 WL's only
    warndlg('Please select only TWO wavelengths to plot!','Error');
end


% --- Executes on button press in OffsetSubtraction.
function OffsetSubtraction_Callback(hObject, eventdata, handles)

interactive = get(handles.InteractiveModeTick,'Value');
handles.SelTraces = [];
switch interactive
    case 1
        timescale = handles.timescale;
        % Select traces interactively, to finish hit RETURN
        handles = SelectTraces(handles);     
        % Check if RAW or CORR and give DATA
        switch handles.rawcorr
            case 'RAW'
                data = transpose(handles.rawsignal);
            case 'CORRECTED'
                data = transpose(handles.corrdata);
        end
        % Get the values to find in the TIME vector
        value = handles.SelTraces(:,2); % column 2 = Time
    case 0
        timescale = handles.timescale;
        handles.SelTraces = inputdlg('Enter space-separated times to plot the transient spectra:',...
             ['Input desired times (in ',timescale,' )'], [1 50]);
        handles.SelTraces = transpose(str2num(handles.SelTraces{:}));
        handles.SelTraces = sort(handles.SelTraces,1); 
        % Check if RAW or CORR and give DATA
        switch handles.rawcorr
            case 'RAW'
                data = transpose(handles.rawsignal);
            case 'CORRECTED'
                data = transpose(handles.corrdata);
        end
         % Get the values to find in the TIME vector
        value = handles.SelTraces(:,1); % column 1 = Time (the only one)
end
% ONCE EVERYTHING IS THERE, THEN DO STUFF
% k = index of the closest match for each selected trace
k = findClosestId2Val(handles.delays,transpose(value));
k = unique(k,'stable');
L = length(k);
% Ask the user whether to define a region for the offset or to use the
% total avg. as a "baseline"
OffsetAvg = questdlg('Do you want to define a baseline or use the spectral average?', ...
    'Offset definition','Average','Define region','Do nothing','Average');
% Calculate the offset-corrected spectra
switch OffsetAvg
    case 'Average'
        data_offset = data - mean(data,1);
    case 'Define region'
        prompt          = {'Enter first and last wavelengths to define offset region:'};
        dlg_title       = 'Define offset region';
        num_lines       = [1 60];
        n_pix           = length(handles.cmprobe);
        defaultans      = {[num2str(handles.cmprobe(n_pix*0.75),4) ' ' num2str(handles.cmprobe(n_pix),'%.3g')]};
        selectoffset    = inputdlg(prompt,dlg_title,num_lines,defaultans);
        OffsetWLs       = str2num(selectoffset{:});
        index           = findClosestId2Val(handles.cmprobe,OffsetWLs);
        % Subtract the offset
        data_offset     = data - mean(data(index(1):index(2),:),1);
end
ydata   = data_offset(:,k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the traces
% Create a new figure with consistent format
% Create a new figure with consistent format
fh = figure();
fh.Position(3)  = 880;
fh.Position(4)  = 425;
fh.Color        = [1 1 1];

% Define the axes
axes2 = axes('Parent',fh);
axes(axes2);

cm=colormap(othercolor('Mrainbow',L));
for n=1:L
   if handles.delays(k(n)) >= 1000
       newdelays = handles.delays(k(n))/1000;
       switch handles.timescale
           case 'fs'
               newtimescale = 'ps';
           case 'ps'
               newtimescale = 'ns';
           case 'ns'
               newtimescale = ['\mu' 's'];
       end
   elseif handles.delays(k(n)) < 1
       newdelays = handles.delays(k(n))*1000;
       switch handles.timescale
           case 'ps'
               newtimescale = 'fs';
           case 'ns'
               newtimescale = 'ps';
       end
   else
       newdelays       = handles.delays(k(n));
       newtimescale    = handles.timescale;
   end
   plot(axes2,handles.cmprobe,ydata(:,n),'LineWidth',1.5,'Marker','none','color',cm(n,:),'DisplayName',[num2str(newdelays,'%.3g') ' ' newtimescale]);
   hold on
end
% Nice formatting
set(axes2,'FontSize',12)
xlabel('Wavenumbers (cm^{-1})','FontWeight','bold','FontSize',13)
ylabel('\DeltaAbs (mOD)','FontWeight','bold','FontSize',13)
axis tight
% Make the legend and hide every second entry if the number of plots is >15
legend('show');
legend('boxoff');
legend('Location','bestoutside');
if L >= 15
    for i=1:2:L
        axes2.Children(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end

% Title and refline
title({handles.datafilename;'OFFSET-CORRECTED DATA'},'Interpreter','none','FontSize',12)
hline = refline(0,0); hline.Color = [0.5 0.5 0.5];
set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
guidata(hObject,handles)

axes2.Units     = 'pixels';
axes2.Position  = [75 60 680 320];
% Tick length
axes2.TickLength    = [0.015 0.035];
% Make the axes resizable
axes2.Units         = 'normalized';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make contour plot of the offset-subtracted data
fh2             = figure();

% Get plot ranges
mintime     = min(handles.delays);
maxtime     = max(handles.delays);
minabs      = min(data_offset(:));
maxabs      = max(data_offset(:));
zminmax     = round(max([abs(minabs) abs(maxabs)]),3);
minwl       = min(handles.cmprobe);
maxwl       = max(handles.cmprobe);
Ncontours   = 40; % 40 contours by default is OK
oldplotranges           = handles.plotranges;
handles.plotranges      = [mintime maxtime minwl maxwl zminmax Ncontours];
% Make the plot
newaxes             = axes('Parent',fh2);
plot2D(handles,transpose(data_offset),newaxes,'Off')
% Restore plot ranges
handles.plotranges      = oldplotranges;
% Formatting
% Make the plot format constant
fh2.Units           = 'pixels';
fh2.Position(2)     = fh2.Position(2).*0.7;
fh2.Position(3)     = 590; % Width in pixels
fh2.Position(4)     = 510; % Height in pixels
fh2.Color           = [1 1 1];

newaxes.FontSize    = 12;
newaxes.Units       = 'pixels';
newaxes.Position    = [75 75 420 385];
% Make uniform, consistent format
newaxes.LineWidth   = 1;
newaxes.TickLength  = [0.015 0.035];

% Make the axes resizable
newaxes.Units       = 'normalized';
% Add correct title
title({handles.datafilename;['DATA w/o OFFSET';'']},'Interpreter','none','FontSize',10)



function percent_whites_Callback(hObject, eventdata, handles)
handles.percentwhites = str2double(handles.percent_whites.String);
switch handles.rawcorr
    case 'RAW'
        plot2D(handles,handles.rawsignal,handles.axes1,'On')
    case 'CORRECTED'
        plot2D(handles,handles.corrdata,handles.axes1,'On')
end
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function percent_whites_CreateFcn(hObject, eventdata, handles)
% hObject    handle to percent_whites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SymmetricColourRange_tick.
function SymmetricColourRange_tick_Callback(hObject, eventdata, handles)
% hObject    handle to SymmetricColourRange_tick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch handles.rawcorr
    case 'RAW'
        plot2D(handles,handles.rawsignal,handles.axes1,'On')
    case 'CORRECTED'
        plot2D(handles,handles.corrdata,handles.axes1,'On')
end
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of SymmetricColourRange_tick


% --- Executes on selection change in SlowMod_selector.
function SlowMod_selector_Callback(hObject, eventdata, handles)

SELECTEDSPECTRUM = handles.SlowMod_selector.Value;
datafilename = handles.datafilename;
% Load the data (replace handles) and update the main handles
        error=0; handles.Nscans = NaN;
        try
           switch handles.DataTypeMenu.Value
               case 1 % Lab 2: UV-Vis Pump - IR Probe (ns)
                   handles = LoadDataIR2(handles);
                   tempdir = strcat(handles.CurrDir.String,filesep,datafilename,filesep,'temp');
                   if exist(tempdir{1},'dir') == 7
                       cd(tempdir{1});
                       filelist = dir('*signal*.csv');
                       handles.Nscans = floor(length(filelist)./2);
                       handles.Nscans_number.String = num2str(handles.Nscans);
                   else
                       handles.Nscans_number.String = 'N/A';
                   end
                   set(handles.ErrorText,'String', '');
               case 2 % Lab 1: UV-Vis/IR Pump - IR probe (fs)
                   handles = LoadDataIR1(handles,SELECTEDSPECTRUM);
                   % TEMP dir story
                   tempdir = strcat(handles.CurrDir.String,filesep,datafilename,filesep,'temp');
                   if exist(tempdir{1},'dir') == 7
                       handles.Nscans_number.String = num2str(handles.Nscans);
                   else
                       handles.Nscans_number.String = 'N/A';
                   end
                   % Clear error string
                   set(handles.ErrorText,'String', '');
           end
        catch err
           set(handles.ErrorText,'String', string(err.message));
           error=1;
        end
        if error == 0
                set(handles.axes1,'Visible','On')
                % Initialise defaults and update the Edits for the first time:
                handles = UpdateEdits(handles);
                guidata(hObject,handles)
                % Override default w/current status of linlog tick
                switch handles.linlogtick.Value
                    case 0
                        handles.linlog = 'lin';
                    case 1
                        handles.linlog = 'log';
                end
        % Plot the data in the main window (for preview)
        % Takes into account the current status of BkgSubTick
        if exist('handles.rawcorr','var') == 0
            switch handles.BkgSubTick.Value
                case 1
                    handles.rawcorr='CORRECTED';
                case 0
                    handles.rawcorr='RAW';
            end
        end
        switch handles.rawcorr
            case 'RAW'
                plot2D(handles,handles.rawsignal,handles.axes1,'On');
            case 'CORRECTED'
                plot2D(handles,handles.corrdata,handles.axes1,'On');
        end
        % Update handles
        guidata(hObject,handles)
        end
% Hints: contents = cellstr(get(hObject,'String')) returns SlowMod_selector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SlowMod_selector


% --- Executes during object creation, after setting all properties.
function SlowMod_selector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SlowMod_selector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GoTodayDir.
function GoTodayDir_Callback(hObject, eventdata, handles)
% Go to the folder
datedir     = datestr(now,'yyyymmdd');
todaydir    = [handles.defaultdir filesep datedir];
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


function plot_skiplevels_Callback(hObject, eventdata, handles)
switch handles.rawcorr
    case 'RAW'
        plot2D(handles,handles.rawsignal,handles.axes1,'On')
    case 'CORRECTED'
        plot2D(handles,handles.corrdata,handles.axes1,'On')
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function plot_skiplevels_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BinScans_Callback(hObject, eventdata, handles)
% Do nothing

% --- Executes during object creation, after setting all properties.
function BinScans_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_colourscheme_Callback(hObject, eventdata, handles)
switch handles.rawcorr
    case 'RAW'
        plot2D(handles,handles.rawsignal,handles.axes1,'On')
    case 'CORRECTED'
        plot2D(handles,handles.corrdata,handles.axes1,'On')
end

function plot_colourscheme_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in AnisotropyCalc.
function AnisotropyCalc_Callback(hObject, eventdata, handles)
% Ask the user what to plot
anisotropy_options = {'UV anisotropy','IR anisotropy','UV magic angle','IR magic angle'};
[anisotropy_typeindx,doplot] = listdlg('ListString',anisotropy_options,'OKstring','Plot','SelectionMode','single','ListSize',[150,120],'PromptString','Select slice type to plot:');

if doplot == 0
    return
end

AnisotropyType  = anisotropy_options{anisotropy_typeindx};
datafilename    = handles.datafilename;
% Load the data (replace handles) and update the main handles
        error=0; handles.Nscans = NaN;
%         try
           switch handles.DataTypeMenu.Value
               case 1 % Lab 2: UV-Vis Pump - IR Probe (ns)
                   handles = LoadDataIR2(handles);
                   tempdir = strcat(handles.CurrDir.String,filesep,datafilename,filesep,'temp');
                   if exist(tempdir{1},'dir') == 7
                       cd(tempdir{1});
                       filelist = dir('*signal*.csv');
                       handles.Nscans = floor(length(filelist)./2);
                       handles.Nscans_number.String = num2str(handles.Nscans);
                   else
                       handles.Nscans_number.String = 'N/A';
                   end
                   set(handles.ErrorText,'String', '');
               case 2 % Lab 1: UV-Vis/IR Pump - IR probe (fs)
                   handles = LoadDataIR1(handles,[],AnisotropyType);
                   % TEMP dir story
                   tempdir = strcat(handles.CurrDir.String,filesep,datafilename,filesep,'temp');
                   if exist(tempdir{1},'dir') == 7
                       handles.Nscans_number.String = num2str(handles.Nscans);
                   else
                       handles.Nscans_number.String = 'N/A';
                   end
                   % Clear error string
                   set(handles.ErrorText,'String', '');
           end
%         catch err
%            set(handles.ErrorText,'String', string(err.message));
%            error=1;
%         end
        if error == 0
                set(handles.axes1,'Visible','On')
                % Initialise defaults and update the Edits for the first time:
                handles = UpdateEdits(handles);
                guidata(hObject,handles)
                % Override default w/current status of linlog tick
                switch handles.linlogtick.Value
                    case 0
                        handles.linlog = 'lin';
                    case 1
                        handles.linlog = 'log';
                end
        % Plot the data in the main window (for preview)
        % Takes into account the current status of BkgSubTick
        if exist('handles.rawcorr','var') == 0
            switch handles.BkgSubTick.Value
                case 1
                    handles.rawcorr='CORRECTED';
                case 0
                    handles.rawcorr='RAW';
            end
        end
        switch handles.rawcorr
            case 'RAW'
                plot2D(handles,handles.rawsignal,handles.axes1,'On');
            case 'CORRECTED'
                plot2D(handles,handles.corrdata,handles.axes1,'On');
        end
        % Update handles
        guidata(hObject,handles)
        end

