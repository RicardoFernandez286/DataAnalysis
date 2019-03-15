function varargout = FitWindow(varargin)
%     Normal usage: 
%     Output = FitWindow(Data_x,Data_y,TraceLabels,YLabel,timescale,FitType)
%               FitType = 1 for normal kinetics, = 2 for SVD (global fit)
%
%     The figure provides a graphical interface to enter the parameters for 
%     up to 4 eponential decays (Tau + amplitude), together with the FWHM and
%     t0 of the response, and an infinite component (step rise) and an offset before t0.
%     The response function is modelled as a normalised Gaussian.
%     
%   Ricardo Fernández-Terán
%   v1.5 / 2018-10-05

% Last Modified by GUIDE v2.5 31-Oct-2017 21:24:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FitWindow_OpeningFcn, ...
                   'gui_OutputFcn',  @FitWindow_OutputFcn, ...
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


% --- Executes just before FitWindow is made visible.
function FitWindow_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FitWindow (see VARARGIN)

% Choose default command line output for FitWindow
handles.output = hObject;

% Disable warning of negative X limits
warning('off','MATLAB:Axes:NegativeLimitsInLogAxis');

% Initialise default values
handles.exp=1;
handles.FitDone=0;
% Linear plot by default
    handles.LinLog=0; 
    handles.LinLog_tick.Value=0;
% Read the traffic lights ;)
handles.OffLight=handles.StatusLight.UserData;
handles.RedLight=handles.DoneButton.UserData;
handles.YellowLight=handles.Preview.UserData;
handles.GreenLight=handles.SaveButton.UserData;

% Read the data (XYY form assumed)
handles.fittype = varargin{6};          % Double [1=Fit traces, 2=Fit SVD]
handles.Data_x = varargin{1};           % 1 x m Array
handles.Data_y = varargin{2};           % m x n Array
handles.YLabel = varargin{4};           % Char
handles.timescale = varargin{5};        % Char

switch handles.fittype
    case 1
        handles.TraceLabels = varargin{3};              % Char (n x 6 for cm-1)
    case 2
        handles.SingVals = varargin{3};                 % Array (containing the numbers of the SingVals
        handles.TraceLabels = {};
        for i=1:length(handles.SingVals)
            handles.TraceLabels{i} = ['SVD Component ' num2str(handles.SingVals(i))];
        end
end

handles.XLabel=['Delays (' handles.timescale ')'];      % Char

% Trace selector initialisation
switch handles.fittype
    case 1
        zcm=[];
        for i=1:size(handles.TraceLabels,1)
            zcm=[zcm; handles.TraceLabels(i,:) ' cm-1'];
        end
        handles.SelectTraces.String = zcm;
        handles.SelectTraces.Value = 1;
        % Plot the incoming data (XYY form assumed, first trace plotted)
        plot(handles.dataAxes,handles.Data_x,handles.Data_y(:,1),'or')
    case 2
        handles.SelectTraces.String = '(SVD global fit)';
        handles.SelectTraces.Value = 1;
        handles.SelectTraces.Enable = 'Off';
        plot(handles.dataAxes,handles.Data_x,handles.Data_y,'o','Linewidth',1.5)
end

% Add labels
xlabel(handles.resAxes,handles.XLabel,'FontSize',13,'FontWeight','bold');
ylabel(handles.dataAxes,handles.YLabel,'FontSize',13,'FontWeight','bold');
ylabel(handles.resAxes,['Res. ' handles.YLabel],'FontSize',13,'FontWeight','bold');
% Formatting
box(handles.dataAxes,'on');
box(handles.resAxes,'on');

% Set limits (tight) and link axes
linkaxes([handles.dataAxes,handles.resAxes],'x');
axis(handles.dataAxes,'tight');
handles.ZL_xlim = [min(handles.Data_x) max(handles.Data_x)];

% Zero lines, new style -(Faster?)
yline(handles.dataAxes,0,'Color',[0.5 0.5 0.5],'HandleVisibility','off');
yline(handles.resAxes,0,'Color',[0.5 0.5 0.5],'HandleVisibility','off');

% Make nice legend @ dataAxes, removing entry for last point (=zeroline)
legend(handles.dataAxes,handles.TraceLabels);
legend(handles.dataAxes,'boxoff'); legend(handles.dataAxes,'Location','northeast')

% Update handles structure
guidata(hObject, handles);

% Wait
uiwait(hObject)

% --- Outputs from this function are returned to the command line.
function varargout = FitWindow_OutputFcn(hObject, eventdata, handles) 

if handles.FitDone==1
    % Output: [FitParams,FitParams_SD,FitFunc,residuals,action]
    varargout{1}=transpose(handles.fitP);
    varargout{2}=handles.fitP_SD;
    varargout{3}=handles.fitfunc;
    varargout{4}=handles.residuals;
    varargout{5}=handles.FitDone;
    varargout{6}=handles.tau_index;
    varargout{7}=handles.c_index;
    varargout{8}=handles.InputParams;
else
    varargout{5}=handles.FitDone;
end
% The figure can be deleted now
delete(handles.FitWindow);



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
handles.FitDone=1;
handles.exit=1;
guidata(hObject,handles);

switch handles.fittype     
    case 1
        % Change the status light
        handles.StatusLight.CData=handles.YellowLight;
        % Build the fit model
        fitparameters = ReadInput(handles); handles.InputParams = fitparameters;
          [handles.fitfunc,handles.p_start,handles.LB,handles.UB,handles.p_text,tau_index,c_index] = BuildFitModel(fitparameters,handles.fittype);
        handles.tau_index=tau_index;
        handles.c_index=c_index;
        % Options for fitting
        options = optimoptions('lsqcurvefit',...
            'MaxFunctionEvaluations',10000,...
            'MaxIterations',2000,...
            'Algorithm','trust-region-reflective',...  %'levenberg-marquardt'
            'OptimalityTolerance',1e-10,...
            'FunctionTolerance',1e-10,...
            'StepTolerance',1e-10);
        % Do the fit
        [handles.fitP,resnorm,handles.residuals,exitflag,output,lambda_fit,jacobian_fit] = lsqcurvefit(handles.fitfunc,handles.p_start,handles.Data_x,handles.Data_y(:,handles.SelectTraces.Value),handles.LB,[],options);
        guidata(hObject,handles);
        fitP=handles.fitP;
        % Update status light
        switch exitflag
            case {1;2;3;4} % The fit worked
                handles.StatusLight.CData=handles.GreenLight;
            case {0;-1;-2} % The fit failed
                handles.StatusLight.CData=handles.RedLight;
        end
    case 2
        % Change the status light
        handles.StatusLight.CData=handles.YellowLight;
        % Build the fit model
        fitparameters = ReadInput(handles);
        handles.InputParams = fitparameters;
        [handles.fitfunc,handles.p_start,handles.LB,handles.UB,handles.p_text,tau_index,c_index] = BuildFitModel(fitparameters,handles.fittype);
        handles.tau_index=tau_index;
        handles.c_index=c_index;
        % Prepare input for NLINMULTIFIT
        x_cell={};
        y_cell={};
        for q=1:length(handles.SingVals)
            x_cell{q} = handles.Data_x;
            y_cell{q} = handles.Data_y(:,q);
        end
            % handles.fitfunc = Already a cell array of function handles
            % handles.p_start = Vector w/ initial guess of fitted params.
        % Fitting Options
        options = statset('MaxIter',2000,...
            'TolFun',1e-10,...
            'TolX',1e-10);
        % Do the fit
        [handles.fitP,RAWresiduals,jacobian_fit,Sigma,MSE] = nlinmultifit(x_cell, y_cell, handles.fitfunc, handles.p_start,options);
        % Change the status light
        handles.StatusLight.CData=handles.GreenLight;
        % Split the residuals into an array (they come out straight)
        handles.residuals=zeros(length(handles.Data_x),length(handles.SingVals));
        for q=1:length(handles.SingVals)
            start=(q-1)*length(handles.Data_x)+1;
            stop=q*length(handles.Data_x);
            handles.residuals(:,q)=RAWresiduals(start:stop);
        end
        % Update status light
% %         switch exitflag
% %             case {1;2;3;4} % The fit worked
% %                 handles.StatusLight.CData=handles.GreenLight;
% %             case {0;-1;-2} % The fit failed
% %                 handles.StatusLight.CData=handles.RedLight;
% %         end
        
        fitP=handles.fitP;
        % Update handles
        guidata(hObject,handles);
        
end

% Update the plot with the fit
handles = UpdatePlot(handles);
guidata(hObject,handles);

if handles.fittype == 1
    %%%%%%%%%%%%%%%%%%%%%%% Calculate statistics
    RES = handles.residuals;
    SSR = sum(RES.^2);
    SSR_s = ['SSR = ' num2str(SSR)];
    CHI2 = sum(RES.^2./abs(handles.Data_y(:,handles.SelectTraces.Value)));
    CHI2_s = [sprintf('\x3c7\x00B2') ' = ' num2str(CHI2)];
    ITER= output.iterations;
    ITER_s = ['Iterations = ' num2str(ITER,'%i')];
    STEP = output.stepsize;
    STEP_s = ['Step size = ' num2str(STEP,2)];
    FUNCEVAL = output.funcCount;
    FUNCEVAL_s =['Fun. evals. = ' num2str(FUNCEVAL,'%i')];
    handles.Statistics_text.String = {SSR_s; '' ; CHI2_s; ''; ITER_s; FUNCEVAL_s; '';STEP_s};
end

% Calculate fit parameter standard deviations
handles.fitP_CI = nlparci(handles.fitP,handles.residuals,'jacobian',jacobian_fit);
fitP_CI=handles.fitP_CI;
fitP_SD=(fitP_CI(:,2)-fitP_CI(:,1))./2;
handles.fitP_SD=fitP_SD;

% Show info about the fit
% The model has been defined before, when doing the fit (fitparameters)
switch handles.fittype
    case 1
        handles = ReadFitOutput(fitparameters,handles,fitP,fitP_SD);
        handles.FitOutput_indicator.String = handles.FitOutput_text;
    case 2
        %%%% TO BE DONE
end

% Update handles
guidata(hObject,handles);

% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles)
% Pass parameters and update handles
handles.FitDone=0;
handles.exit=1;
guidata(hObject,handles);

% The figure can be deleted now
uiresume

% --- Executes on selection change in NumExp_menu.
function NumExp_menu_Callback(hObject, eventdata, handles)
contents    = cellstr(get(hObject,'String'));
exps        = contents{get(hObject,'Value')};
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
        eval(['handles.Tau' num2str(i) '_value.Enable = ' off ';']);
        eval(['handles.c' num2str(i) '_value.Enable = ' off ';']);
        eval(['handles.Tau' num2str(i) '_hold.Enable = ' off ';']);
        eval(['handles.c' num2str(i) '_hold.Enable = ' off ';']);
        eval(['handles.Beta' num2str(i) '_value.Enable = ' off ';']);
        eval(['handles.Beta' num2str(i) '_hold.Enable = ' off ';']);
end
% Switch on the selected Nº of exponentials
if handles.exp ~= 0
    for i=1:handles.exp
        eval(['handles.Tau' num2str(i) '_value.Enable = ' on ';']);
        eval(['handles.c' num2str(i) '_value.Enable = ' on ';']);
        eval(['handles.Tau' num2str(i) '_hold.Enable = ' on ';']);
        eval(['handles.c' num2str(i) '_hold.Enable = ' on ';']);
        eval(['handles.Beta' num2str(i) '_value.Enable = ' on ';']);
        eval(['handles.Beta' num2str(i) '_hold.Enable = ' on ';']);
    end
end
% Update handles
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function NumExp_menu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Preview.
function Preview_Callback(hObject, eventdata, handles)
% Pass the parameters
handles.FitDone=2;
handles.exit=1;

handles=UpdatePlot(handles);

% Update handles
guidata(hObject,handles);

% --- Executes on button press in ResetZoom.
function ResetZoom_Callback(hObject, eventdata, handles)
xlim(handles.dataAxes,[min(handles.Data_x) max(handles.Data_x)]);
set(handles.dataAxes,'XScale','lin');
set(handles.resAxes,'XScale','lin');
handles.ZL_xlim = [min(handles.Data_x) max(handles.Data_x)];

% Zero lines, new style -(Faster?)
yline(handles.dataAxes,0,'Color',[0.5 0.5 0.5],'HandleVisibility','off');
yline(handles.resAxes,0,'Color',[0.5 0.5 0.5],'HandleVisibility','off');

% Update handles
guidata(hObject,handles);

% --- Executes on button press in LinLog_tick.
function LinLog_tick_Callback(hObject, eventdata, handles)
handles.LinLog = get(hObject,'Value');
switch handles.LinLog
    case 1 % Log
        set(handles.resAxes,'XScale','log');
        set(handles.dataAxes,'XScale','log');
        handles.ZL_xlim = [min(handles.Data_x(handles.Data_x>0)) max(handles.Data_x(handles.Data_x>0))];
    case 0 % Lin
        set(handles.resAxes,'XScale','lin');
        set(handles.dataAxes,'XScale','lin');
        handles.ZL_xlim = [min(handles.Data_x) max(handles.Data_x)];
end

% Zero lines, new style -(Faster?)
yline(handles.dataAxes,0,'Color',[0.5 0.5 0.5],'HandleVisibility','off');
yline(handles.resAxes,0,'Color',[0.5 0.5 0.5],'HandleVisibility','off');

% Update handles
guidata(hObject,handles);

% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Create new figure and axes
fh = figure();
dataAxes2 = axes('Parent',fh,'Position',[64, 218, 721, 405],'Units','normalized');
resAxes2 = axes('Position',[64, 64, 721, 122]);
switch handles.FitDone
    case 1 % Real fit done
        p_values = handles.fitP;
        % Get the values of the function and plot it
        syms t;
        previewplot = handles.fitfunc(p_values,t);
        axes(dataAxes2);
        limits = handles.ZL_xlim;
        % Redo the plots
        plot(handles.Data_x,handles.Data_y(:,handles.SelectTraces.Value),'or')
        hold on;
        fplot(dataAxes2,previewplot,limits,'b','LineWidth',2);
        hold off;
        % Plot the residuals from the fit
        plot(resAxes,handles.Data_x,handles.residuals,'or');
        xlim(dataAxes2,limits);
    case 2 % Preview
        % Read the inputs
        Parameters = ReadInput(handles);
        % Build the model equation
        [fitfunc,p_start,LB,p_text] = BuildFitModel(Parameters,handles.fittype);
        % Plot the model equation with the given parameters (initial guess)
        syms t;
        p_values = p_start;
        previewplot = fitfunc(p_values,t);
        axes(dataAxes2);
        limits = handles.ZL_xlim;
        % Redo the plots
        plot(handles.Data_x,handles.Data_y(:,handles.SelectTraces.Value),'or')
        hold on;
        fplot(dataAxes2,previewplot,limits,'b','LineWidth',2);
        hold off;
        xlim(dataAxes2,limits);
        % Plot the residuals from the initial guess
        plot(resAxes2,handles.Data_x,handles.Data_y(:,handles.SelectTraces.Value) - fitfunc(p_start,handles.Data_x),'or');
        xlim(resAxes2,limits);
end

% Add labels
xlabel(resAxes2,handles.XLabel,'FontSize',14,'FontWeight','bold');
ylabel(dataAxes2,handles.YLabel,'FontSize',14,'FontWeight','bold');
ylabel(resAxes2,['Res. ' handles.YLabel],'FontSize',14,'FontWeight','bold');

% Lin/Log
switch handles.LinLog
    case 1
        set(resAxes2,'XScale','log');
        set(dataAxes2,'XScale','log');
    case 0
        set(resAxes2,'XScale','lin');
        set(dataAxes2,'XScale','lin');
end

% Zero lines, new style -(Faster?)
yline(handles.dataAxes,0,'Color',[0.5 0.5 0.5],'HandleVisibility','off');
yline(handles.resAxes,0,'Color',[0.5 0.5 0.5],'HandleVisibility','off');

%%% INTERNAL FUNCTIONS
function Output = ReadInput(handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The output has the following structure:
%     Tau_values   | Tau1 Tau2 Tau3 Tau4
%     Tau_holds    | Tau1 Tau2 Tau3 Tau4
%     C_values     | c1   c2   c3   c4   cInf  Offset
%     C_holds      | c1   c2   c3   c4   cInf  Offset
%     Param_values | FWHM t0
%     Param_holds  | FWHM t0
%     Switches     | GaussIRF  StretchedExp (1 or 0 = On/Off)
%     FitDone        | 0=Cancel, 1=Fit, 2=Preview
%     NumExp       | Number of exponentials to fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocate memory
Tau_values  = zeros(1,4);
Tau_holds   = zeros(1,4);
Beta_values  = zeros(1,4);
Beta_holds   = zeros(1,4);
switch handles.fittype
    case 1
        C_values    = zeros(1,6);
        C_holds     = zeros(1,6);
    case 2
        n=length(handles.SingVals);
        C_values    = zeros(n,6);
        C_holds     = zeros(n,6);    
end
Param_values= zeros(1,2);
Param_holds = zeros(1,2);
% Give output
if handles.FitDone ~= 0 %User did not press "Cancel" button
% Copy special parameters
    % Read
    FWHM=str2num(handles.FWHM_value.String);
    FWHM_h=handles.FWHM_hold.Value;
    t0=str2num(handles.t0_value.String);
    t0_h=handles.t0_hold.Value;
    Offset = str2num(handles.Offset_value.String);
    Offset_h = handles.Offset_hold.Value;
    % Write
    Param_values = [FWHM t0];
    Param_holds  = [FWHM_h t0_h];
    C_values(:,6) = transpose(Offset);
    C_holds(:,6)= transpose(Offset_h);
% Get the infinite
    % Read
    cInf=str2num(handles.Infty_value.String);
    cInf_h=handles.Infty_hold.Value;
    % Write
    C_values(:,5)=transpose(cInf);
    C_holds(:,5)=transpose(cInf_h);
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
        currBeta=['Beta' num2str(i) '_value'];
        currC       =   ['c' num2str(i) '_value'];
        currTauhold =   ['Tau' num2str(i) '_hold'];
        currBetahold=['Beta' num2str(i) '_hold'];
        currChold   =   ['c' num2str(i) '_hold'];
        Tau_values(i)  =   str2num(eval(['handles.' currTau '.String']));
        Tau_holds(i)   =   eval(['handles.' currTauhold '.Value']);
        Beta_values(i)  =   str2num(eval(['handles.' currBeta '.String']));
        Beta_holds(i)   =   eval(['handles.' currBetahold '.Value']);
        C_values(:,i)    =   transpose(str2num(eval(['handles.' currC '.String'])));
        C_holds(:,i)     =   transpose(eval(['handles.' currChold '.Value']));
    end
end
    Output = struct;
    Output.Tau_values = Tau_values;
    Output.Tau_holds = Tau_holds;
    Output.Beta_values = Beta_values;
    Output.Beta_holds = Beta_holds;
    Output.C_values= C_values;
    Output.C_holds= C_holds;
    Output.Param_values = Param_values;
    Output.Param_holds = Param_holds;
    Output.FitDone = handles.FitDone;
    Output.NumExp = handles.exp;
    Output.Switches = [handles.GaussIRF_tick.Value handles.StretchedExp_tick.Value];
    
function handles = UpdatePlot(handles)
n_plotpoints = 10000;
switch handles.FitDone
    case 1 % REAL FIT DONE
        % Read the inputs
        Parameters = ReadInput(handles);
        % Build the model equation
        [fitfunc,~,~,~] = BuildFitModel(Parameters,handles.fittype);
        % Plot the model equation with the given parameters (initial guess)
        syms t;
        p_values = handles.fitP;
        if handles.fittype==2
            nFits=length(handles.SingVals);
            cmap=colormap(lines(nFits));
        else
            nFits=1;
        end
        cla(handles.dataAxes);
        cla(handles.resAxes);
        for FitNo=1:nFits
            if handles.fittype==2
                plotfunc = fitfunc{FitNo};
            else
                plotfunc = fitfunc;
            end
            limits = handles.ZL_xlim;
            % Redo the plots
            hold(handles.dataAxes,'on');
            hold(handles.resAxes,'on');
            switch handles.LinLog
                case 1
                    timeaxis  = logspace(log10(limits(1)),log10(limits(2)),n_plotpoints);
                case 0
                    timeaxis  = linspace(limits(1),limits(2),n_plotpoints);
            end
            fitted_data   = double(plotfunc(p_values,timeaxis));
            if handles.fittype==2
                dp(FitNo) = plot(handles.dataAxes,handles.Data_x,handles.Data_y(:,FitNo),'o','color',cmap(FitNo,:),'LineWidth',0.025);
                fp(FitNo) = plot(handles.dataAxes,timeaxis,fitted_data,'-','LineWidth',2,'color',cmap(FitNo,:),'DisplayName',['SVD Component ' num2str(FitNo)]);
                rp(FitNo) = plot(handles.resAxes,handles.Data_x,(handles.Data_y(:,FitNo) - plotfunc(p_values,handles.Data_x)),'-o','color',cmap(FitNo,:));
            else
                plot(handles.dataAxes,handles.Data_x,handles.Data_y(:,handles.SelectTraces.Value),'o','color',[0.2 0.2 1])               
                plot(handles.dataAxes,timeaxis,fitted_data,'-','LineWidth',2,'color','b')
                plot(handles.resAxes,handles.Data_x,handles.Data_y(:,handles.SelectTraces.Value) - plotfunc(p_values,handles.Data_x),'-o','color',[0.2 0.2 1]);
            end
            xlim(handles.dataAxes,limits);
            xlim(handles.resAxes,limits);
        end
    case 2 % Preview
        % Read the inputs
        Parameters = ReadInput(handles);
        % Build the model equation
        [fitfunc,p_start,LB,p_text] = BuildFitModel(Parameters,handles.fittype);
        % Plot the model equation with the given parameters (initial guess)
        syms t;
        p_values = p_start;
        if handles.fittype==2
            nFits=length(handles.SingVals);
            cmap=colormap(lines(nFits));
        else
            nFits=1;
        end
        cla(handles.dataAxes);
        cla(handles.resAxes);
        for FitNo=1:nFits
            if handles.fittype==2
                plotfunc = fitfunc{FitNo};
            else
                plotfunc = fitfunc;
            end
            previewplot = plotfunc(p_values,t);
            limits = handles.ZL_xlim;
            % Redo the plots
            hold(handles.dataAxes,'on');
            hold(handles.resAxes,'on')
            switch handles.LinLog
                case 1
                    timeaxis  = logspace(log10(limits(1)),log10(limits(2)),n_plotpoints);
                case 0
                    timeaxis  = linspace(limits(1),limits(2),n_plotpoints);
            end
            fitted_data   = double(plotfunc(p_values,timeaxis));
            if handles.fittype==2
                dp(FitNo) = plot(handles.dataAxes,handles.Data_x,handles.Data_y(:,FitNo),'o','color',cmap(FitNo,:),'LineWidth',0.025);
                fp(FitNo) = plot(handles.dataAxes,timeaxis,fitted_data,'-','LineWidth',2,'color',cmap(FitNo,:),'DisplayName',['SVD Component ' num2str(FitNo)]);
                rp(FitNo) = plot(handles.resAxes,handles.Data_x,(handles.Data_y(:,FitNo) - plotfunc(p_values,handles.Data_x)),'-o','color',cmap(FitNo,:));
            else
                plot(handles.dataAxes,handles.Data_x,handles.Data_y(:,handles.SelectTraces.Value),'o','color',[0.2 0.2 1])
                plot(handles.dataAxes,timeaxis,fitted_data,'-','LineWidth',2,'color','b')
                plot(handles.resAxes,handles.Data_x,handles.Data_y(:,handles.SelectTraces.Value) - plotfunc(p_values,handles.Data_x),'-o','color',[0.2 0.2 1]);
            end
            xlim(handles.dataAxes,limits);
            xlim(handles.resAxes,limits);
        end
    case 3
        % REDO THE PLOT WHEN THE WAVELENGTH IS CHANGED
        axes(handles.dataAxes);
        limits = handles.ZL_xlim;
        plot(handles.Data_x,handles.Data_y(:,handles.SelectTraces.Value),'-o','color',[1 0.2 0.2])
        xlim(handles.dataAxes,limits);
        cla(handles.resAxes);
end

% Add labels
xlabel(handles.resAxes,handles.XLabel,'FontSize',14,'FontWeight','bold');
ylabel(handles.dataAxes,handles.YLabel,'FontSize',14,'FontWeight','bold');
ylabel(handles.resAxes,['Res. ' handles.YLabel],'FontSize',14,'FontWeight','bold');

% Lin/Log
switch handles.LinLog
    case 1
        set(handles.resAxes,'XScale','log');
        set(handles.dataAxes,'XScale','log');
    case 0
        set(handles.resAxes,'XScale','lin');
        set(handles.dataAxes,'XScale','lin');
end

% Zero lines, new style -(Faster?)
yline(handles.dataAxes,0,'Color',[0.5 0.5 0.5],'HandleVisibility','off');
yline(handles.resAxes,0,'Color',[0.5 0.5 0.5],'HandleVisibility','off');


if handles.fittype == 2
    % Make nice legend @ dataAxes, removing entry for last point (=zeroline)
    legend(handles.dataAxes,'boxoff'); legend(handles.dataAxes,'Location','northeast')
    legend(fp(1:nFits))
else
    legend(handles.dataAxes,'off');
end

% --- Executes on selection change in SelectTraces.
function SelectTraces_Callback(hObject, eventdata, handles)
if handles.fittype ~= 2
    % Update the plot with the fit
    guidata(hObject,handles);
    handles.FitDone= 3;
    handles = UpdatePlot(handles);
end

% --- Executes during object creation, after setting all properties.
function SelectTraces_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectTraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in DoneButton.
function DoneButton_Callback(hObject, eventdata, handles)
% Get output and close the figure
uiresume

% --- Executes on button press in StatusLight.
function StatusLight_Callback(hObject, eventdata, handles)
%Do nothing


% --------------------------------------------------------------------
function ZoomIn_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ZoomIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zoom
zoom(gcf,'Direction','in')


% --------------------------------------------------------------------
function ZoomOut_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ZoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zoom
zoom(gcf,'Direction','out')

% --------------------------------------------------------------------
function Pan_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to Pan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pan


% --------------------------------------------------------------------
function DataCursor_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to DataCursor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datacursormode


% --------------------------------------------------------------------
function InsertLegend_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to InsertLegend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function Beta1_value_Callback(hObject, eventdata, handles)
% hObject    handle to Beta1_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Beta1_value as text
%        str2double(get(hObject,'String')) returns contents of Beta1_value as a double


% --- Executes during object creation, after setting all properties.
function Beta1_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Beta1_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Beta1_hold.
function Beta1_hold_Callback(hObject, eventdata, handles)
% hObject    handle to Beta1_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Beta1_hold



function Beta2_value_Callback(hObject, eventdata, handles)
% hObject    handle to Beta2_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Beta2_value as text
%        str2double(get(hObject,'String')) returns contents of Beta2_value as a double


% --- Executes during object creation, after setting all properties.
function Beta2_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Beta2_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Beta2_hold.
function Beta2_hold_Callback(hObject, eventdata, handles)
% hObject    handle to Beta2_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Beta2_hold



function Beta3_value_Callback(hObject, eventdata, handles)
% hObject    handle to Beta3_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Beta3_value as text
%        str2double(get(hObject,'String')) returns contents of Beta3_value as a double


% --- Executes during object creation, after setting all properties.
function Beta3_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Beta3_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Beta3_hold.
function Beta3_hold_Callback(hObject, eventdata, handles)
% hObject    handle to Beta3_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Beta3_hold



function Beta4_value_Callback(hObject, eventdata, handles)
% hObject    handle to Beta4_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Beta4_value as text
%        str2double(get(hObject,'String')) returns contents of Beta4_value as a double


% --- Executes during object creation, after setting all properties.
function Beta4_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Beta4_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Beta4_hold.
function Beta4_hold_Callback(hObject, eventdata, handles)
% hObject    handle to Beta4_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Beta4_hold


% --- Executes on button press in GaussIRF_tick.
function GaussIRF_tick_Callback(hObject, eventdata, handles)
% hObject    handle to GaussIRF_tick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tick=get(hObject,'Value');
if tick == 1 % GaussIRF implies NO beta
    opt='on'; opt2='off';
    handles.StretchedExp_tick.Value = 0;
    handles.Beta_text.Visible = 'off';
    handles.HoldBeta_text.Visible = 'off';
    for i=1:4
        eval(['handles.Beta' num2str(i) '_value.Visible = opt2;']);
        eval(['handles.Beta' num2str(i) '_hold.Visible = opt2;']);
    end
elseif tick == 0 % No GaussIRF, don't touch anything
    opt='off';
end
% Switch GaussIRF items on/off
handles.FWHM_text.Visible = opt;
handles.FWHM_value.Enable = opt;
handles.FWHM_hold.Enable = opt;
handles.FWHM_value.Visible = opt;
handles.FWHM_hold.Visible = opt;

% Switch the betas
if handles.StretchedExp_tick.Value == 1
    opt='on';
elseif handles.StretchedExp_tick.Value == 0
    opt='off';
end

% --- Executes on button press in StretchedExp_tick.
function StretchedExp_tick_Callback(hObject, eventdata, handles)
% hObject    handle to StretchedExp_tick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tick=get(hObject,'Value');
if tick == 1
    opt='on';
    % Switch the IRF off
    handles.GaussIRF_tick.Value = 0;
    handles.FWHM_text.Visible = 'off';
    handles.FWHM_value.Enable = 'off';
    handles.FWHM_hold.Enable = 'off';
    handles.FWHM_value.Visible = 'off';
    handles.FWHM_hold.Visible = 'off';
elseif tick == 0 % No beta, don't touch IRF
    opt='off';
end

% Switch the Betas on or off
handles.Beta_text.Visible = opt;
handles.HoldBeta_text.Visible = opt;
for i=1:4
    eval(['handles.Beta' num2str(i) '_value.Visible = opt;']);
    eval(['handles.Beta' num2str(i) '_hold.Visible = opt;']);
end
