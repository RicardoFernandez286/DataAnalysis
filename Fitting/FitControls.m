function varargout = FitControls(varargin)
% FitControls: MATLAB code for FitControls.fig
%     The figure provides a graphical interface to enter the parameters for 
%     up to 4 eponential decays (Tau + amplitude), together with the FWHM and
%     t0 of the response, and an infinite component (step rise) and an offset before t0.
%     The response function is modelled as a normalised Gaussian.
%     
%     Normal usage: Output = FitControls()
%    
%     Output (given upon pressing either "Fit!" or "Cancel" is a structure
%     with the following fields:
%         Tau_values      Values of the lifetimes
%         Tau_holds       Switches to hold (or not) the lifetimes
%         C_values        Values of the amplitudes. C_values(5) = Infinite
%         C_holds         Switches to hold (or not) the amplitudes
%         Param_values    Contains the extra params: FWHM t0 Offset (in that order)
%         Param_holds     Switches to hold (or not) the corresp. parameters
%         DoFit           Switch which tells whether the user pressed "Fit" or "Cancel"
%
%   Ricardo Fernández-Terán
%   v1.0 / 2017-06-07

% Last Modified by GUIDE v2.5 15-Jun-2017 16:28:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FitControls_OpeningFcn, ...
                   'gui_OutputFcn',  @FitControls_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before FitControls is made visible.
function FitControls_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FitControls (see VARARGIN)

% Choose default command line output for FitControls
handles.output = hObject;

% Initialise default fields
handles.exp=1;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FitControls wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = FitControls_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The output has the following structure:
% varargout{i}
%     1               Tau1 Tau2 Tau3 Tau4
%     3   HOLDS for   Tau1 Tau2 Tau3 Tau4
%     2               c1   c2   c3   c4   cInf
%     4   HOLDS for   c1   c2   c3   c4   cInf
%     5               FWHM t0   Offset
%     6   HOLDS for   FWHM t0   Offset
%     7               doFit (1 or 0)
%     8               NumExp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocate memory
Tau_values  = zeros(1,4);
Tau_holds   = zeros(1,4);
C_values    = zeros(1,5);
C_holds     = zeros(1,5);
Param_values= zeros(1,3);
Param_holds = zeros(1,3);
% Give output
if handles.dofit ~= 0 %User did not press "Cancel" button
% Copy special parameters
    % Read
    FWHM=str2num(handles.FWHM_value.String);
    FWHM_h=handles.FWHM_hold.Value;
    t0=str2num(handles.t0_value.String);
    t0_h=handles.t0_hold.Value;
    Offset = str2num(handles.Offset_value.String);
    Offset_h = handles.Offset_hold.Value;
    % Write
    Param_values = [FWHM t0 Offset];
    Param_holds  = [FWHM_h t0_h Offset_h];
% Get the infinite
    % Read
    cInf=str2num(handles.Infty_value.String);
    cInf_h=handles.Infty_hold.Value;
    % Write
    C_values(5)=cInf;
    C_holds(5)=cInf_h;
% Check number of exponentials
    exps = handles.NumExp_menu.String{handles.NumExp_menu.Value};
    switch exps
        case 'One exponential'
            handles.exp=1;
        case 'Two exponentials'
            handles.exp=2;
        case 'Three exponentials'
            handles.exp=3;
        case 'Four exponentials'
            handles.exp=4;
        case 'Step rise'
            handles.exp=0;
    end
    for i=1:handles.exp
        currTau=['Tau' num2str(i) '_value'];
        currC       =   ['c' num2str(i) '_value'];
        currTauhold =   ['Tau' num2str(i) '_hold'];
        currChold   =   ['c' num2str(i) '_hold'];
        Tau_values(i)  =   str2num(eval(['handles.' currTau '.String']));
        C_values(i)    =   str2num(eval(['handles.' currC '.String']));
        Tau_holds(i)   =   eval(['handles.' currTauhold '.Value']);
        C_holds(i)     =   eval(['handles.' currChold '.Value']);
    end
end
    Output = struct;
    Output.Tau_values = Tau_values;
    Output.Tau_holds = Tau_holds;
    Output.C_values= C_values;
    Output.C_holds= C_holds;
    Output.Param_values = Param_values;
    Output.Param_holds = Param_holds;
    Output.DoFit = handles.dofit;
    Output.NumExp = handles.exp;
    varargout{1}= Output;
% The figure can be deleted now
delete(handles.figure1);



function FWHM_value_Callback(hObject, eventdata, handles)
% hObject    handle to FWHM_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FWHM_value as text
%        str2double(get(hObject,'String')) returns contents of FWHM_value as a double


% --- Executes during object creation, after setting all properties.
function FWHM_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FWHM_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FWHM_hold.
function FWHM_hold_Callback(hObject, eventdata, handles)
% hObject    handle to FWHM_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FWHM_hold



function t0_value_Callback(hObject, eventdata, handles)
% hObject    handle to t0_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t0_value as text
%        str2double(get(hObject,'String')) returns contents of t0_value as a double


% --- Executes during object creation, after setting all properties.
function t0_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t0_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in t0_hold.
function t0_hold_Callback(hObject, eventdata, handles)
% hObject    handle to t0_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of t0_hold



function Tau1_value_Callback(hObject, eventdata, handles)
% hObject    handle to Tau1_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tau1_value as text
%        str2double(get(hObject,'String')) returns contents of Tau1_value as a double


% --- Executes during object creation, after setting all properties.
function Tau1_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tau1_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Tau1_hold.
function Tau1_hold_Callback(hObject, eventdata, handles)
% hObject    handle to Tau1_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Tau1_hold



function Tau2_value_Callback(hObject, eventdata, handles)
% hObject    handle to Tau2_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tau2_value as text
%        str2double(get(hObject,'String')) returns contents of Tau2_value as a double


% --- Executes during object creation, after setting all properties.
function Tau2_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tau2_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Tau2_hold.
function Tau2_hold_Callback(hObject, eventdata, handles)
% hObject    handle to Tau2_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Tau2_hold



function Tau3_value_Callback(hObject, eventdata, handles)
% hObject    handle to Tau3_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tau3_value as text
%        str2double(get(hObject,'String')) returns contents of Tau3_value as a double


% --- Executes during object creation, after setting all properties.
function Tau3_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tau3_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Tau3_hold.
function Tau3_hold_Callback(hObject, eventdata, handles)
% hObject    handle to Tau3_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Tau3_hold



function Tau4_value_Callback(hObject, eventdata, handles)
% hObject    handle to Tau4_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tau4_value as text
%        str2double(get(hObject,'String')) returns contents of Tau4_value as a double


% --- Executes during object creation, after setting all properties.
function Tau4_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tau4_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Tau4_hold.
function Tau4_hold_Callback(hObject, eventdata, handles)
% hObject    handle to Tau4_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Tau4_hold



function c1_value_Callback(hObject, eventdata, handles)
% hObject    handle to c1_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c1_value as text
%        str2double(get(hObject,'String')) returns contents of c1_value as a double


% --- Executes during object creation, after setting all properties.
function c1_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c1_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in c1_hold.
function c1_hold_Callback(hObject, eventdata, handles)
% hObject    handle to c1_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c1_hold



function c2_value_Callback(hObject, eventdata, handles)
% hObject    handle to c2_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c2_value as text
%        str2double(get(hObject,'String')) returns contents of c2_value as a double


% --- Executes during object creation, after setting all properties.
function c2_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c2_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in c2_hold.
function c2_hold_Callback(hObject, eventdata, handles)
% hObject    handle to c2_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c2_hold



function c3_value_Callback(hObject, eventdata, handles)
% hObject    handle to c3_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c3_value as text
%        str2double(get(hObject,'String')) returns contents of c3_value as a double


% --- Executes during object creation, after setting all properties.
function c3_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c3_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in c3_hold.
function c3_hold_Callback(hObject, eventdata, handles)
% hObject    handle to c3_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c3_hold



function c4_value_Callback(hObject, eventdata, handles)
% hObject    handle to c4_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c4_value as text
%        str2double(get(hObject,'String')) returns contents of c4_value as a double


% --- Executes during object creation, after setting all properties.
function c4_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c4_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in c4_hold.
function c4_hold_Callback(hObject, eventdata, handles)
% hObject    handle to c4_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c4_hold



function Infty_value_Callback(hObject, eventdata, handles)
% hObject    handle to Infty_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Infty_value as text
%        str2double(get(hObject,'String')) returns contents of Infty_value as a double


% --- Executes during object creation, after setting all properties.
function Infty_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Infty_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Infty_hold.
function Infty_hold_Callback(hObject, eventdata, handles)
% hObject    handle to Infty_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Infty_hold



function Offset_value_Callback(hObject, eventdata, handles)
% hObject    handle to Offset_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Offset_value as text
%        str2double(get(hObject,'String')) returns contents of Offset_value as a double


% --- Executes during object creation, after setting all properties.
function Offset_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Offset_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Offset_hold.
function Offset_hold_Callback(hObject, eventdata, handles)
% hObject    handle to Offset_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Offset_hold


% --- Executes on button press in DoFit.
function DoFit_Callback(hObject, eventdata, handles)
% hObject    handle to DoFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dofit=1;
handles.exit=1;
handles.outstring='Do Fit!';
guidata(hObject,handles);
% The GUI is no longer waiting, just close it via opening function
uiresume(handles.figure1)

% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dofit=0;
handles.exit=1;
handles.outstring='Cancel';
guidata(hObject,handles);
% The GUI is no longer waiting, just close it via opening function
uiresume(handles.figure1)

% --- Executes on selection change in NumExp_menu.
function NumExp_menu_Callback(hObject, eventdata, handles)
% hObject    handle to NumExp_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
exps = contents{get(hObject,'Value')};
switch exps
    case 'One exponential'
        handles.exp=1;
    case 'Two exponentials'
        handles.exp=2;
    case 'Three exponentials'
        handles.exp=3;
    case 'Four exponentials'
        handles.exp=4;
    case 'Step rise'
        handles.exp=0;
end
on=char('on'); off=char('off');
% Switch everything off
for i=1:4
        currTau=['Tau' num2str(i) '_value'];
        currC=['c' num2str(i) '_value'];
        currTauhold=['Tau' num2str(i) '_hold'];
        currChold=['c' num2str(i) '_hold'];
        eval(['handles.' currTau '.Enable = ' off ';']);
        eval(['handles.' currC '.Enable = ' off ';']);
        eval(['handles.' currTauhold '.Enable = ' off ';']);
        eval(['handles.' currChold '.Enable = ' off ';']);
end
% Switch on the selected Nº of exponentials
if handles.exp ~= 0
    for i=1:handles.exp
        currTau=['Tau' num2str(i) '_value'];
        currC=['c' num2str(i) '_value'];
        currTauhold=['Tau' num2str(i) '_hold'];
        currChold=['c' num2str(i) '_hold'];
        eval(['handles.' currTau '.Enable = ' on ';']);
        eval(['handles.' currC '.Enable = ' on ';']);
        eval(['handles.' currTauhold '.Enable = ' on ';']);
        eval(['handles.' currChold '.Enable = ' on ';']);
    end
end
% Update handles
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function NumExp_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumExp_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Preview.
function Preview_Callback(hObject, eventdata, handles)
% hObject    handle to Preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dofit=2;
handles.exit=1;
handles.outstring='Preview';
guidata(hObject,handles);
% The GUI is no longer waiting, just close it via opening function
uiresume(handles.figure1)