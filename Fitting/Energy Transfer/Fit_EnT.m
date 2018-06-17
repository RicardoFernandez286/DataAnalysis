function [Taus,Amps,Tau_err,Amp_err,SSR] = Fit_EnT(time,FW_data,BW_data,start_params,LB,UB,ShowOutput,DoFit,varargin)
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
% Ricardo Fernández-Terán, v0.9a - 2018.06.04

%% Build the fit model
% Set the options for fitting
options = optimoptions('lsqcurvefit',...
            'MaxFunctionEvaluations',10000,...
            'MaxIterations',10000,...
            'Algorithm','trust-region-reflective',...  %'levenberg-marquardt' 'trust-region-reflective'
            'OptimalityTolerance',1e-5,...
            'FunctionTolerance',1e-5,...
            'StepTolerance',1e-5);

% Define the parameter locations
Tau_index   = [1 2 3 4];
Amp_index   = [5 6 7 8]; 

% Get the input parameters
% start_params(Tau_index) = start_params(Tau_index);
% start_params(Amp_index) = start_amps;

%% Parse the input data
% Split the arrays
diagFW      = FW_data(:,1);
diagBW      = BW_data(:,1);
xpeakFW     = FW_data(:,2);
xpeakBW     = BW_data(:,2);

if DoFit==1
%% Fit the data
    [fitted_param,resnorm,residuals,exitflag,output,lambda,jacobian_fit] = lsqcurvefit(@EnT_kinetics,start_params,time,[diagFW,xpeakFW,diagBW,xpeakBW],LB,UB,options);
    plot_param      = fitted_param;
    fitted_param    = fitted_param';

%% Calculate parameters and parameter errors
    fitted_CI   = nlparci(fitted_param,residuals,'jacobian',jacobian_fit);
    fitted_CI   = (fitted_CI(:,2) - fitted_CI(:,1))/2;

    Taus        = fitted_param(Tau_index);
    Amps        = fitted_param(Amp_index);
    Tau_err     = fitted_CI(Tau_index);
    Amp_err     = fitted_CI(Amp_index);

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
            '\n' ...
            'AMPLITUDES: \n' ...
            '       A_0        = %.4g ' char(177) ' %.4g \n' ...
            '       AB_0       = %.4g ' char(177) ' %.4g \n' ...
            '       B_0        = %.4g ' char(177) ' %.4g \n' ...
            '       BA_0       = %.4g ' char(177) ' %.4g \n' ...
            '\n'];             
%         Output = [Taus Tau_err]';
        Output = [[Taus Tau_err];[Amps Amp_err]]';   
        fprintf(formatspec,Output(:));
    end
else
    plot_param = start_params;
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
plot(ax_FW,time,10.*xpeakFW,'^r','LineWidth',1,'DisplayName','Re(^{13}CO) \rightarrow Re(^{12}CO) \rm{\times10}');
plot(ax_BW,time,10.*xpeakBW,'vb','LineWidth',1,'DisplayName','Re(^{13}CO) \leftarrow Re(^{12}CO) \rm{\times10}');

% Set axes limits
axis(ax_FW,'tight');
axis(ax_BW,'tight');

%% Plot the fit results
% Plot the fit/initial guess kinetics (evaluating the function)
t_values    = linspace(time(1),time(end),500);
fitted_data = EnT_kinetics(plot_param,t_values);

fit_diagFW  = fitted_data(:,1);
fit_xpeakFW = 10.*fitted_data(:,2);
fit_diagBW  = fitted_data(:,3);
fit_xpeakBW = 10.*fitted_data(:,4);

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
    text(ax_FW,time(end),0.6,['\tau_{ET,fw} = ' num2str(Taus(3),'%.3g') ' \pm ' num2str(Tau_err(3),'%.3g') ' ps'],'HorizontalAlignment','right','FontSize',12,'FontWeight','bold');
    text(ax_BW,time(end),0.6,['\tau_{ET,bw} = ' num2str(Taus(4),'%.3g') ' \pm ' num2str(Tau_err(4),'%.3g') ' ps'],'HorizontalAlignment','right','FontSize',12,'FontWeight','bold');
end

% Make figure resizable
fh.Units        = 'normalized';