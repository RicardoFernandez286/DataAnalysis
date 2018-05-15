%close all
% minmax = [-5 300];
% syms t
% %parameters = FitControls;
% DoFit = parameters.DoFit;
% if DoFit ==1
% Tau_values = parameters.Tau_values;
% Tau_holds = parameters.Tau_holds;
% C_values = parameters.C_values;
% C_holds = parameters.C_holds;
% Param_values = parameters.Param_values;
% Param_holds = parameters.Param_holds;
% else
% end
% 
% FWHM = Param_values(1);
% sigma=FWHM/(2*sqrt(2*log(2)));
% t0   = Param_values(2);

% y0 = Param_values(3);
% y1 = C_values(1)*0.5*exp(sigma^2/(2*Tau_values(1)^2)).*exp(-(t-t0)/Tau_values(1)).*(1+erf((t-t0-sigma^2/Tau_values(1))/sigma/sqrt(2)));
% y2 = C_values(2)*0.5*exp(sigma^2/(2*Tau_values(2)^2)).*exp(-(t-t0)/Tau_values(2)).*(1+erf((t-t0-sigma^2/Tau_values(2))/sigma/sqrt(2)));
% y3 = C_values(3)*0.5*exp(sigma^2/(2*Tau_values(3)^2)).*exp(-(t-t0)/Tau_values(3)).*(1+erf((t-t0-sigma^2/Tau_values(3))/sigma/sqrt(2)));
% y4 = C_values(4)*0.5*exp(sigma^2/(2*Tau_values(4)^2)).*exp(-(t-t0)/Tau_values(4)).*(1+erf((t-t0-sigma^2/Tau_values(4))/sigma/sqrt(2)));
% yInf = C_values(5)*0.5*(1+erf((t-t0)/sigma/sqrt(2)));
% Y = y0 + y1 + y2 + y3 + y4 + yInf;
% 
% 
% fplot(Y,minmax,'b','LineWidth',2);
% hline = refline(0,0);  hline.Color = [0.5 0.5 0.5];
% 
% noise=0.02;
% [x,y]=fplot(Y,minmax);
% hold on
% ynoise=y+noise*randn(size(y));
% plot(x,ynoise,'-o')
% hold off

% syms cOff c1 c2 c3 c4 cInf tau1 tau2 tau3 tau4 fwhm t0 
% syms c tau
% 
% sigma=fwhm/(2*sqrt(2*log(2)));
% 
% %f0 = cOff;
% f1 = c(1)*0.5*exp(sigma^2/(2*tau(1)^2)).*exp(-(t-t0)/tau(1)).*(1+erf((t-t0-sigma^2/tau(1))/sigma/sqrt(2)));
% f2 = c(2)*0.5*exp(sigma^2/(2*tau(2)^2)).*exp(-(t-t0)/tau(2)).*(1+erf((t-t0-sigma^2/tau(2))/sigma/sqrt(2)));
% %f3 = c(3)*0.5*exp()^2/(2*tau(3)^2)).*exp(-(t-t0)/tau(3)).*(1+erf((t-t0-(FWHM/(2*sqrt(2*log(2))))^2/tau(3))/(FWHM/(2*sqrt(2*log(2))))/sqrt(2)));
% %f4 = c(4)*0.5*exp(sigma^2/(2*tau(4)^2)).*exp(-(t-t0)/tau(4)).*(1+erf((t-t0-sigma^2/tau(4))/sigma/sqrt(2)));
% %fInf = cInf*0.5*(1+erf((t-t0)/sigma/sqrt(2)));
% 
% FUN= f1 + f2;
% 
% lsqcurvefit()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gather input
% Query parameters
parameters = FitControls;
% Read given parameters
DoFit = parameters.DoFit;

Tau_values = parameters.Tau_values;
Tau_holds = parameters.Tau_holds;
C_values = parameters.C_values;
C_holds = parameters.C_holds;
Param_values = parameters.Param_values;
Param_holds = parameters.Param_holds;
NumExp = parameters.NumExp;
%% Build the model
% Var is a number corresponding to the variable passed later to the fit
% model
var=0;

% p_start is a vector containing the starting values of the parameters
p_start=[];

% p_text is a cell array containing all parameters' names in order
p_text={};

% Check whether FWHM, t0 and Offset are held, then set them as a fixed
% number or as a variable to fit
if Param_holds(1) == 0 % FWHM
    var         = var+1;
    p_i         = num2str(var);
    p_start(var)= Param_values(1);
    FWHM_f      = ['p(' p_i ')'];
    p_text{var} = 'FWHM';
else
    FWHM_f = num2str(Param_values(1));
end

if Param_holds(2) == 0 % t0
    var         = var+1;
    p_i         = num2str(var);
    p_start(var)= Param_values(2);
    t0_f        = ['p(' p_i ')'];
    p_text{var} = 't0';
else
    t0_f = num2str(Param_values(2));
end

if Param_holds(3) == 0 % Offset
    var           = var+1;
    p_i           = num2str(var);
    p_start(var)  = Param_values(3);
    Offset_f      = ['p(' p_i ')'];
    p_text{var}   = 'Offset';
else
    Offset_f = num2str(Param_values(3));
end

% Put in the infinite component
if C_holds(5) == 0 % Infinite component
    var         = var+1;
    p_i         = num2str(var);
    p_start(var)= C_values(5);
    cInf_f      = ['p(' p_i ')'];
    p_text{var} = 'cInf';
else
    cInf_f = num2str(C_values(5));
end
    infmodel =[Offset_f '+' cInf_f '*0.5*(1+erf((t-' t0_f ')/(' FWHM_f '/(2*sqrt(log(2))))))'];

% Build the exponentials
maxExp=4; % Set the maximum number of exponentials - for future upgrades
expmodel=[];
if NumExp ~= 0
    % Preallocate memory and initialise variables
    C = {};
    Tau = {};
    % Check each exponential and write the model
    for i=1:NumExp
        if C_holds(i) == 0 % Amplitude
            var         = var+1;
            p_i         = num2str(var);
            p_start(var)= C_values(i);
            C{i}        = ['p(' p_i ')'];
            p_text{var} = ['C' num2str(i)];
        else
            C{i} = num2str(C_values(i));
        end
        if Tau_holds(i) == 0 % Lifetime
            var         = var+1;
            p_i         = num2str(var);
            p_start(var)= Tau_values(i);
            Tau{i}      = ['p(' p_i ')'];
            p_text{var} = [' tau' num2str(i)];
        else
            Tau{i} = num2str(Tau_values(i));
        end
        % Write the model
        expmodel=strcat(expmodel,'+0.5*',C{i},'.*exp((',FWHM_f,'/(2*sqrt(2*log(2)))^2)/(2*(',Tau{i},')^2)).*exp(-(t-',t0_f,')/',Tau{i},').*(1+erf((t-',t0_f,'-(',FWHM_f,'/(2*sqrt(2*log(2)))^2)/',Tau{i},')/(sqrt(2)*',FWHM_f,'/(2*sqrt(2*log(2))))))');
    end
end
fitfunc_model=strcat(infmodel,expmodel);
fitfunc = eval(['@(p,t) ' fitfunc_model]);

%% Do the fitting
options = optimoptions('lsqcurvefit',...
            'MaxFunctionEvaluations',5000,...
            'MaxIterations',1000,...
            'Algorithm','levenberg-marquardt',...
            'OptimalityTolerance',1e-10,...
            'FunctionTolerance',1e-10,...
            'StepTolerance',1e-10,...
            'UseParallel',true);

[pfit,~,residuals,exitflag,output] = lsqcurvefit(fitfunc,p_start,x,ynoisy,[],[],options)
%% DO STUFF AFTER FIT
switch DoFit
    case 1
        hold on
        plot(x,fitfunc(pfit,x),'b','LineWidth',2)
        plot(x,ynoisy,'-or')
        hold off
    case 2
        syms t;
        guessplot = fitfunc(p_start,t);
        plotrange=[-5 50];

        fplot(guessplot,plotrange,'LineWidth',2)
        [x,y]=fplot(guessplot,plotrange);
        noise=0.05;
        ynoisy=y+noise*randn(size(y));
        hold on
        plot(x,ynoisy,'-or')
        hline = refline(0,0);  hline.Color = [0.5 0.5 0.5];
end
