function varargout = Gaussian2D_setparams(varargin)
%GAUSSIAN2D_SETPARAMS MATLAB code file for Gaussian2D_setparams.fig
%      GAUSSIAN2D_SETPARAMS, by itself, creates a new GAUSSIAN2D_SETPARAMS or raises the existing
%      singleton*.
%
%      H = GAUSSIAN2D_SETPARAMS returns the handle to a new GAUSSIAN2D_SETPARAMS or the handle to
%      the existing singleton*.
%
%      GAUSSIAN2D_SETPARAMS('Property','Value',...) creates a new GAUSSIAN2D_SETPARAMS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to Gaussian2D_setparams_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GAUSSIAN2D_SETPARAMS('CALLBACK') and GAUSSIAN2D_SETPARAMS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GAUSSIAN2D_SETPARAMS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Gaussian2D_setparams

% Last Modified by GUIDE v2.5 18-Mar-2019 18:01:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Gaussian2D_setparams_OpeningFcn, ...
                   'gui_OutputFcn',  @Gaussian2D_setparams_OutputFcn, ...
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


function Gaussian2D_setparams_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for Gaussian2D_setparams
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Gaussian2D_setparams wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function varargout = Gaussian2D_setparams_OutputFcn(hObject, eventdata, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;


function ExitButton_Callback(hObject, eventdata, handles)


function DataTypeMenu_Callback(hObject, eventdata, handles)


function DataTypeMenu_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Population_delay_Callback(hObject, eventdata, handles)


function Population_delay_CreateFcn(hObject, eventdata, handles)

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SaveProbeCal_Callback(hObject, eventdata, handles)


function Transient2D_mode_Callback(hObject, eventdata, handles)


function Transient2D_mode_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function PowerDependence_Callback(hObject, eventdata, handles)


function SubtractSpectra_Callback(hObject, eventdata, handles)


function Fit_Callback(hObject, eventdata, handles)


function ShiftT2_Callback(hObject, eventdata, handles)


function SpectralDiffusion_Callback(hObject, eventdata, handles)


function IntegralDynamics_Callback(hObject, eventdata, handles)


function ShowTimeDomain_Callback(hObject, eventdata, handles)


function ShowPhasing_Callback(hObject, eventdata, handles)


function ShowProbeCalibration_Callback(hObject, eventdata, handles)


function PixelNumber_Callback(hObject, eventdata, handles)


function PixelNumber_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function apodise_method_Callback(hObject, eventdata, handles)


function apodise_method_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function zeropad_tick_Callback(hObject, eventdata, handles)


function zeropad_next2k_Callback(hObject, eventdata, handles)


function zeropad_factor_Callback(hObject, eventdata, handles)


function zeropad_factor_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function phase_Npoints_Callback(hObject, eventdata, handles)


function phase_Npoints_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function phase_method_Callback(hObject, eventdata, handles)


function phase_method_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function phase_coeffs_Callback(hObject, eventdata, handles)


function phase_coeffs_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function make2Dplot_Callback(hObject, eventdata, handles)


function make3Dplot_Callback(hObject, eventdata, handles)


function plot_Kinetics_Callback(hObject, eventdata, handles)


function plot_Slices_Callback(hObject, eventdata, handles)


function InteractiveModeTick_Callback(hObject, eventdata, handles)


function Normalise_Callback(hObject, eventdata, handles)


function plot_Ncontours_Callback(hObject, eventdata, handles)


function plot_Ncontours_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editMaxZ_Callback(hObject, eventdata, handles)


function editMaxZ_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function timescalemenu_Callback(hObject, eventdata, handles)


function timescalemenu_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function plot_Nwhites_Callback(hObject, eventdata, handles)


function plot_Nwhites_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ShowContoursTick_Callback(hObject, eventdata, handles)


function plot_contourLineStyle_Callback(hObject, eventdata, handles)


function plot_contourLineStyle_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function plot_pumpdirection_Callback(hObject, eventdata, handles)


function plot_pumpdirection_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function plot_axislegend_Callback(hObject, eventdata, handles)


function plot_axislegend_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function plot_colorrange_Callback(hObject, eventdata, handles)


function plot_colorrange_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function plot_skiplevels_Callback(hObject, eventdata, handles)


function plot_skiplevels_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function plot_contourfill_Callback(hObject, eventdata, handles)


function plot_colourscheme_Callback(hObject, eventdata, handles)


function plot_colourscheme_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SymColRange_tick_Callback(hObject, eventdata, handles)


function CutPlot_tick_Callback(hObject, eventdata, handles)


function EditProbeAxis_tick_Callback(hObject, eventdata, handles)


function editWLmin_Callback(hObject, eventdata, handles)


function editWLmin_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editWLmax_Callback(hObject, eventdata, handles)


function editWLmax_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ShowSpecDiff_Callback(hObject, eventdata, handles)


function DoFit_Callback(hObject, eventdata, handles)


function Preview_Callback(hObject, eventdata, handles)


function FitOK_Callback(hObject, eventdata, handles)


function SaveButton_Callback(hObject, eventdata, handles)


function AddPeak_Callback(hObject, eventdata, handles)


function RemovePeak_Callback(hObject, eventdata, handles)
