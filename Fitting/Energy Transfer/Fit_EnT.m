function [fitted_param,varargout] = Fit_EnT(time,FW_data,BW_data,start_params,LB,UB,ShowOutput,DoFit,varargin)
% Fitting function used to fit forward<->backward energy transfer problems
% Input:
%    time = vector with the time steps in ps
%    fw   = two-column matrix with the peak intensities for the diagonal and the FW cross-peak
%    bw   = two-column matrix with the peak intensities for the diagonal and the BW cross-peak
%    
%    In this context,
%     FW energy transfer means Re(13CO) to Re(12CO) (i.e. uphill)
%     BW energy transfer means Re(12CO) to Re(13CO) (i.e. downhill)
%
% Ricardo Fern�ndez-Ter�n, v0.9a - 2018.06.04
fit_routine = 'Analytical LS'; % 'Numerical LS' 'Analytical LS' 'Analytical GA'

%% Build the fit model
% Set the options for fitting
options = optimoptions('lsqcurvefit',...
            'MaxFunctionEvaluations',1000,...
            'MaxIterations',1000,...
            'Algorithm','trust-region-reflective',...  %'levenberg-marquardt' 'trust-region-reflective'
            'OptimalityTolerance',1e-15,...
            'FunctionTolerance',1e-15,...
            'StepTolerance',1e-15);%,...
%             'PlotFcn',{@optimplotresnorm,@optimplotfirstorderopt});

% Define the parameter locations
Tau_index   = 1:4;
Amp_index   = 5; 

% Get the input parameters
% start_params(Tau_index) = start_params(Tau_index);
% start_params(Amp_index) = start_amps;

%% Parse the input data
% Split the arrays
diagFW      = FW_data(:,1);
diagBW      = BW_data(:,1);
xpeakFW     = 10.*FW_data(:,2);
xpeakBW     = 10.*BW_data(:,2);

%% Solve the differential equation analytically and prepare the fit function
kB      = 8.6173303*10^-5; % eV K^-1
Temp    = 298; % K
dW      = 46; % cm-1
ev2cm   = 8065.54;
factor  = exp(-dW./((kB*Temp*ev2cm)));

[A_func,aB_func,bA_func,B_func] = EnT_kinetics_analytic(time,FW_data,BW_data);

%% Define the fitting and evaluation functions

function valuesFW = eval_fitfuncFW(P,t)
    A_val   = double(A_func(P(1),P(2),P(3),P(4),P(5),t));
    aB_val  = double(aB_func(P(1),P(2),P(3),P(4),P(5),t));
    valuesFW= [A_val,10.*aB_val];
end

function valuesBW = eval_fitfuncBW(P,t)
    B_val   = double(B_func(P(1),P(2),P(3),P(4),P(5),t));
    bA_val  = double(bA_func(P(1),P(2),P(3),P(4),P(5),t));
    valuesBW= [B_val,10.*bA_val];
end

function valuesALL = eval_fitfuncALL(P,t)
    A_val   = double(A_func(P(1),P(2),P(3),P(4),P(5),t));
    aB_val  = double(aB_func(P(1),P(2),P(3),P(4),P(5),t));
    B_val   = double(B_func(P(1),P(2),P(3),P(4),P(5),t));
    bA_val  = double(bA_func(P(1),P(2),P(3),P(4),P(5),t));
    valuesALL = [A_val,aB_val,B_val,bA_val];
end

%% Define the global optimisation routines (GA, particle swarm, etc.)
function fitted_param = runGAfit()
    rng default
%     options = optimoptions('ga','UseParallel', true, 'UseVectorized', false);
    fitted_param = ga(@obj_fun,5,[],[],[],[],LB,UB,[]);
    % Nested function to calculate the objectivefunciton
    function ssr = obj_fun(P)
        ssr = sumsqr(eval_fitfuncALL(P,time)-[diagFW,xpeakFW,diagBW,xpeakBW]);
    end
end

%% Fit the data
if DoFit==1
tic

switch fit_routine
    case 'Numerical LS'
% Numerical solution of the ODE system
    [fitted_param,resnorm,residuals,exitflag,output,lambda,jacobian_fit] = lsqcurvefit(@EnT_kinetics,start_params,time,[diagFW,xpeakFW,xpeakBW,diagBW],LB,UB,options);

    case 'Analytical LS'
% Analytical solution of the ODE system
    [fitted_param,resnorm,residuals,exitflag,output,lambda,jacobian_fit] = lsqcurvefit(@eval_fitfuncALL,start_params,time,[diagFW,xpeakFW,diagBW,xpeakBW],LB,UB,options);

    case 'Analytical GA'
% Genetic algorithm fitting using the analytical solution of the ODE system
    fitted_param = runGAfit;
end

toc
%% Calculate parameters and parameter errors
plot_param      = fitted_param;
fitted_param    = fitted_param';

    Taus        = fitted_param(Tau_index);
    Amps        = fitted_param(Amp_index);
%    Amps        = NaN(length(Amps),1);

switch fit_routine
    case 'Analytical GA'
        Tau_err     = NaN(length(Taus),1);
        Amp_err     = NaN(length(Amps));
    otherwise
        fitted_CI   = nlparci(fitted_param,residuals,'jacobian',jacobian_fit);
        fitted_CI   = (fitted_CI(:,2) - fitted_CI(:,1))/2;
        Tau_err     = fitted_CI(Tau_index);
        Amp_err     = fitted_CI(Amp_index);
end

%% Give fit results
    if ShowOutput == 1
        formatspec = [...
            '\n' ...
            'FIT RESULTS: \n' ...
            '=============== \n' ...
            '\n' ...
            'LIFETIMES: \n' ...
            '       Tau_fast   = %.4g ' char(177) ' %.4g ps \n' ...
            '       Tau_slow   = %.4g ' char(177) ' %.4g ps \n' ...
            '       Tau_FW     = %.4g ' char(177) ' %.4g ps \n' ...
            '       Tau_BW     = %.4g ' char(177) ' %.4g ps \n' ...
            'AMPLITUDES: \n' ...
            '       C          = %.4g ' char(177) ' %.4g \n' ...
            '\n' ...
            ];             
        Output = [[Taus,Tau_err];[Amps,Amp_err]]';
%         Output = [[Taus];[Amps]]';   

        fprintf(formatspec,Output(:));
    end
else
    plot_param = start_params;
end


varargout{1} = eval_fitfuncALL(plot_param,time);

if ShowOutput == 0
    return
end
%% Plot the data
% Create figure
if isempty(varargin) == 0 && isvalid(varargin{1}) == 1
    fh              = varargin{1};
    clf(fh,'reset');
else
    fh              = figure;
    fh.Units        = 'normalized';
    fh.Position(2)  = 0.1;
    fh.Position(4)  = fh.Position(4)*2;
end

fh.Color        = [1 1 1];
fh.Units        = 'pixels';
ax_FW           = subplot(2,1,1);
ax_BW           = subplot(2,1,2);
ax_FW.FontSize  = 12;
ax_BW.FontSize  = 12;

box(ax_FW,'on');
box(ax_BW,'on');

hold(ax_FW,'on');
hold(ax_BW,'on');

% Plot the diagonal peaks
plot(ax_FW,time,diagFW,'or','LineWidth',1,'DisplayName','Diagonal Re(^{13}CO)');
plot(ax_BW,time,diagBW,'ob','LineWidth',1,'DisplayName','Diagonal Re(^{12}CO)');

% Plot the cross peaks
plot(ax_FW,time,xpeakFW,'^r','LineWidth',1,'DisplayName','Re(^{13}CO) \rightarrow Re(^{12}CO) \rm{\times10}');
plot(ax_BW,time,xpeakBW,'vb','LineWidth',1,'DisplayName','Re(^{13}CO) \leftarrow Re(^{12}CO) \rm{\times10}');

% Set axes limits
axis(ax_FW,'tight');
axis(ax_BW,'tight');

%% Plot the fit results
% Plot the fit/initial guess kinetics (evaluating the function)
t_values    = (linspace(time(1),time(end),500))';
% fitted_data = EnT_kinetics(plot_param,t_values);
fitted_data = eval_fitfuncALL(plot_param,t_values);

fit_diagFW  = fitted_data(:,1);
fit_xpeakFW = fitted_data(:,2);

fit_diagBW  = fitted_data(:,3);
fit_xpeakBW = fitted_data(:,4);

plot(ax_FW,t_values,fit_diagFW,'-','LineWidth',1.5,'color','r','DisplayName','Diagonal - FW (Fit)');
plot(ax_FW,t_values,fit_xpeakFW,'-','LineWidth',1.5,'color','r','DisplayName','Re(^{13}CO) \rightarrow Re(^{12}CO) (fit)');
plot(ax_BW,t_values,fit_diagBW,'-','LineWidth',1.5,'color','b','DisplayName','Diagonal - BW (Fit)');
plot(ax_BW,t_values,fit_xpeakBW,'-','LineWidth',1.5,'color','b','DisplayName','Re(^{13}CO) \leftarrow Re(^{12}CO) (fit)');

hold(ax_FW,'off');
hold(ax_BW,'off');

%% Set plot formatting
% Titles
title(ax_FW,'Forward energy transfer','FontSize',13);
title(ax_BW,'Backward energy transfer','FontSize',13);

% Axis labels
xlabel(ax_FW,'t_2 delay (ps)','FontWeight','bold','FontSize',12);
xlabel(ax_BW,'t_2 delay (ps)','FontWeight','bold','FontSize',12);

ylabel(ax_FW,'Intensity (a.u.)','FontWeight','bold','FontSize',12);
ylabel(ax_BW,'Intensity (a.u.)','FontWeight','bold','FontSize',12);

% Axis limits
xlim(ax_FW,[-5 time(end)+5]);
xlim(ax_BW,[-5 time(end)+5]);

ylim(ax_FW,[-0.1 1.1]);
ylim(ax_BW,[-0.1 1.1]);

% Legends
lh_FW = legend(ax_FW,'show');
legend(ax_FW,'boxoff','FontWeight','bold')
legend(ax_FW,{},'FontWeight','bold')

lh_BW = legend(ax_BW,'show');
legend(ax_BW,'boxoff')
legend(ax_BW,{},'FontWeight','bold')

% Delete legends of the fits
lh_FW.String(3:4) = [];
lh_BW.String(3:4) = [];

% Add annotations
if DoFit==1
    text(ax_FW,time(end),0.6,['\tau_{EnT,fw} = ' num2str(Taus(3),'%.3g') ' \pm ' num2str(Tau_err(3),'%.3g') ' ps'],'HorizontalAlignment','right','FontSize',12,'FontWeight','bold');
    text(ax_BW,time(end),0.6,['\tau_{EnT,bw} = ' num2str(Taus(4),'%.3g') ' \pm ' num2str(Tau_err(4),'%.3g') ' ps'],'HorizontalAlignment','right','FontSize',12,'FontWeight','bold');
end

% Add reflines
% hlineFW = refline(ax_FW,[0 0]);
% hlineFW.Color = [0.5 0.5 0.5];
% set(get(get(hlineFW,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
% 
% hlineBW = refline(ax_BW,[0 0]);
% hlineBW.Color = [0.5 0.5 0.5];
% set(get(get(hlineBW,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

% Make figure resizable
fh.Units        = 'normalized';
end