function varargout = Gaussian2D_fitparam(varargin)
% GAUSSIAN2D_FITPARAM MATLAB code for Gaussian2D_fitparam.fig
%      GAUSSIAN2D_FITPARAM, by itself, creates a new GAUSSIAN2D_FITPARAM or raises the existing
%      singleton*.
%
%      H = GAUSSIAN2D_FITPARAM returns the handle to a new GAUSSIAN2D_FITPARAM or the handle to
%      the existing singleton*.
%
%      GAUSSIAN2D_FITPARAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GAUSSIAN2D_FITPARAM.M with the given input arguments.
%
%      GAUSSIAN2D_FITPARAM('Property','Value',...) creates a new GAUSSIAN2D_FITPARAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Gaussian2D_fitparam_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Gaussian2D_fitparam_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Gaussian2D_fitparam

% Last Modified by GUIDE v2.5 21-Aug-2019 20:03:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Gaussian2D_fitparam_OpeningFcn, ...
                   'gui_OutputFcn',  @Gaussian2D_fitparam_OutputFcn, ...
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


% --- Executes just before Gaussian2D_fitparam is made visible.
function Gaussian2D_fitparam_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Gaussian2D_fitparam (see VARARGIN)

% Choose default command line output for Gaussian2D_fitparam
handles.output = hObject;

% Initialise variables
handles.Return = 1;
handles.t2_startfit_text.String = num2str(0.1);

if isempty(varargin)
    handles.FitParam_table.Data =    ...
       [{'1980'}    {'2060'}    {'1980' }    {'2060'  }
        {'25'  }    {'25'  }    {'-2060'}    {'-1980' }
        {'15'  }    {'15'  }    {'-1'   }    {'-1'    }
        {'10'  }    {'10'  }    {'-1'   }    {'-1'    }
        {'1'   }    {'1'   }    {'0h'   }    {'0h'    }
        {'Diag'}    {'Diag'}    {'Xpeak'}    {'Xpeak'}];
    handles.t2_endfit_text.String = num2str(60);
    handles.equal_SxSy.Value        = 1;
    handles.diffSyfor12.Value       = 0;
elseif strcmp(varargin{1},'Re1213 VET')
    handles.FitParam_table.Data =    ...
       [{'1979'}    {'2028'}    {'1979' }    {'2028'  }
        {'11'  }    {'11'  }    {'-2028'}    {'-1979' }
        {'10'  }    {'10'  }    {'-1'   }    {'-1'    }
        {'10'  }    {'10'  }    {'-1'   }    {'-1'    }
        {'1'   }    {'1'   }    {'0h'   }    {'0h'    }
        {'Diag'}    {'Diag'}    {'Xpeak'}    {'Xpeak'}];
    handles.t2_startfit_text.String = num2str(0);
    handles.t2_endfit_text.String   = num2str(60);
    handles.equal_SxSy.Value        = 0;
    handles.diffSyfor12.Value       = 1;
else
    handles.FitParam_table.Data =    ...
   [{'1980'}    {'2060'}    {'1980' }    {'2060'  }
    {'25'  }    {'25'  }    {'-2060'}    {'-1980' }
    {'15'  }    {'15'  }    {'-1'   }    {'-1'    }
    {'10'  }    {'10'  }    {'-1'   }    {'-1'    }
    {'1'   }    {'1'   }    {'0h'   }    {'0h'    }
    {'Diag'}    {'Diag'}    {'Xpeak'}    {'Xpeak'}];
    t2delays = varargin{2};
    handles.t2_startfit_text.String = num2str(min(t2delays(t2delays>0)));
    handles.t2_endfit_text.String   = num2str(max(t2delays));
    handles.equal_SxSy.Value        = 1;
    handles.diffSyfor12.Value       = 0;
end

handles.FitParam_table.ColumnEditable = true(1,size(handles.FitParam_table.Data,2));

% Update handles structure
guidata(hObject, handles);

% Wait
uiwait(handles.Gaussian2D_fitparam)


% --- Outputs from this function are returned to the command line.
function varargout = Gaussian2D_fitparam_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isempty([handles.FitParam_table.Data{:}])
    varargout{1} = {};
elseif handles.Return == 0
    Data    = handles.FitParam_table.Data;
    t2range = [str2double(handles.t2_startfit_text.String), str2double(handles.t2_endfit_text.String)];
    varargout{1} = Data;
    varargout{2} = t2range;
else
    varargout{1} = {};
    varargout{2} = [];
end
% Give default output values
varargout{3} = handles.equal_SxSy.Value;    % equal_SxSy
varargout{4} = handles.diffSyfor12.Value;   % different Sy for the 1-2 transitions

varargout{5} = handles.Return;
close Gaussian2D_fitparam

function AddPeak_Callback(hObject, eventdata, handles)
handles.FitParam_table.Data = [handles.FitParam_table.Data repmat({''},6,1,size(handles.FitParam_table.Data,3))];
handles.FitParam_table.ColumnEditable = true(1,size(handles.FitParam_table.Data,2));
if size(handles.FitParam_table.Data,2) >= 1
    handles.RemovePeak.Enable = 'On';
end
guidata(hObject, handles);

function RemovePeak_Callback(hObject, eventdata, handles)
handles.FitParam_table.Data(:,end,:) = [];
handles.FitParam_table.ColumnEditable = true(1,size(handles.FitParam_table.Data,2));
if size(handles.FitParam_table.Data,2) < 1
    handles.RemovePeak.Enable = 'Off';
end
guidata(hObject, handles);

function DoneButton_Callback(hObject, eventdata, handles)
handles.Return = 0;
guidata(hObject, handles);
uiresume

function CancelButton_Callback(hObject, eventdata, handles)
handles.FitParam_table.Data = {};
handles.Return = 1;
guidata(hObject, handles);
uiresume


function edit2_Callback(hObject, eventdata, handles)


function edit2_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit3_Callback(hObject, eventdata, handles)


function edit3_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function t2_startfit_text_Callback(hObject, eventdata, handles)


function t2_startfit_text_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function t2_endfit_text_Callback(hObject, eventdata, handles)


function t2_endfit_text_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function equal_SxSy_Callback(hObject, eventdata, handles)


function diffSyfor12_Callback(hObject, eventdata, handles)
