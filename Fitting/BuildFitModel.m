function [fitfunc,p_start,LB,UB,p_text,Tau_index,C_index] = BuildFitModel(parameters,fittype)
% This function will build the fit function used to fit the data
% Usage: [fitfunc,p_start,LB,p_text] = BuildFitModel(parameters)
% The input is a structure with the following fields:
%         Tau_values    Tau1 Tau2 Tau3 Tau4
%         Tau_holds     Tau1 Tau2 Tau3 Tau4
%         C_values      c1   c2   c3   c4  cInf Offset
%         C_holds       c1   c2   c3   c4  cInf Offset
%         Param_values  FWHM t0
%         Param_holds   FWHM t0
%         DoFit         0=Cancel, 1=Fit, 2=Preview
%         NumExp        Number of exponentials to fit
%         Switches      GaussIRF  StretchedExp (1 or 0)    
%
% The equation to fit is the result of an offset (c0) added to the 
% convolution of a sum of NumExp and an infinite component
%
%   Ricardo Fernández-Terán
%   v1.2 / 2017-07-24

%% Read the input structure
% Read given parameters
DoFit = parameters.FitDone;

Tau_values = parameters.Tau_values;
Tau_holds = parameters.Tau_holds;

Beta_values = parameters.Beta_values;
Beta_holds = parameters.Beta_holds;

C_values = parameters.C_values;
C_holds = parameters.C_holds;

Param_values = parameters.Param_values;
Param_holds = parameters.Param_holds;

% Read the switches
GaussIRF = parameters.Switches(1);
Stretched = parameters.Switches(2);

% Get the number of exponentials to fit
NumExp = parameters.NumExp;

% Get the number of traces to be simultaneously fitted
if fittype==1
    nFits=1;
elseif fittype==2
    nFits=size(C_values,1);
end
%% Build the model
% Var is a number corresponding to the variable passed later to the fit
% model
var=0;

% p_start is a vector containing the starting values of the parameters
p_start=[];

% p_text is a cell array containing all parameters' names in order
p_text={};

% LB and UB are arrays containing the lower and upper limits for all fitted
% parameters. 
% All time constants must be positive, while the others can go to -Inf
% Beta can only take values between 0 and 1.
LB=[]; UB=[];

% Tau_index and C_index are two arrays containing the index numbers of
% those parameters in the var vector. If held, index=0.
Tau_index=zeros(1,NumExp);
C_index=zeros(nFits,6);

% Check whether FWHM and t0 are held, then set them as a fixed
% number or as a variable to fit
if GaussIRF ~= 0 % IF GAUSSIAN IRF INCLUDED, THEN DO IT
    if Param_holds(1) == 0 % FWHM
        var         = var+1;
        p_i         = num2str(var);
        p_start(var)= Param_values(1);
        FWHM_f      = ['p(' p_i ')'];
        p_text{var} = 'FWHM';
        LB(var) = 0; % FWHM can't be negative!
        UB(var) = +Inf;
    else
        FWHM_f = num2str(Param_values(1));
    end
end

if Param_holds(2) == 0 % t0
    var         = var+1;
    p_i         = num2str(var);
    p_start(var)= Param_values(2);
    t0_f        = ['p(' p_i ')'];
    p_text{var} = 't0';
    LB(var) = -Inf;
    UB(var) = +Inf;
else
    t0_f = num2str(Param_values(2));
end

expmodel=cell(nFits,1);
expmodel_temp=cell(nFits,1);
infmodel=cell(nFits,1);
for FitNo=1:nFits
    %%% Offset and infinite component
    % Put in the Offset (before t0)
    if C_holds(FitNo,6) == 0
        var           = var+1;
        C_index(FitNo,6) = var;
        p_i           = num2str(var);
        p_start(var)  = C_values(FitNo,6);
        Offset_f      = ['p(' p_i ')'];
        p_text{var}   = ['Offset' '.' num2str(FitNo)];
        LB(var) = -Inf;
        UB(var) = +Inf;
    else
        Offset_f = num2str(C_values(FitNo,6));
        C_index(FitNo,6) = 0;
    end

    % Put in the infinite component
    if C_holds(FitNo,5) == 0 % Infinite component
        var         = var+1;
        C_index(FitNo,5) = var;
        p_i         = num2str(var);
        p_start(var)= C_values(FitNo,5);
        cInf_f      = ['p(' p_i ')'];
        p_text{var} = ['cInf' '.' num2str(FitNo)];
        LB(var) = -Inf;
        UB(var) = +Inf;
    else
        cInf_f = num2str(C_values(FitNo,5));
        C_index(FitNo,5) = 0;
    end

    % Write the infinite component equation
    if GaussIRF ~= 0 % GAUSSIAN IRF
        if C_holds(FitNo,5)==1 && C_values(FitNo,5)==0
            infmodel{FitNo}=[];
        else
            infmodel{FitNo}=[Offset_f '+' cInf_f '*0.5*(1+erf((t-' t0_f ')/(' FWHM_f '/(2*sqrt(log(2))))))'];
        end
    else % STEP FUNCTION IRF
        if C_holds(FitNo,5)==1 && C_values(FitNo,5)==0
            infmodel{FitNo}=[];
        else
            infmodel{FitNo}=[Offset_f '+' cInf_f '*heaviside(t-' t0_f ')'];
        end
    end
    
    %%% Build the exponentials
    maxExp=4; % Set the maximum number of exponentials - for future upgrades;

    if NumExp ~= 0
        % Preallocate memory and initialise variables
        C = cell(nFits,NumExp);
        Tau = cell(1,NumExp);
        Beta = cell(1,NumExp);
        % Check each exponential and write the model
        for ExpNo=1:NumExp
            if C_holds(FitNo,ExpNo) == 0 % Amplitude
                var         = var+1;
                C_index(FitNo,ExpNo) = var;
                p_i         = num2str(var);
                p_start(var)= C_values(FitNo,ExpNo);
                C{FitNo,ExpNo}  = ['p(' p_i ')'];
                p_text{var} = ['C' num2str(ExpNo) '.' num2str(FitNo)];
                LB(var) = -Inf;
                UB(var) = +Inf;
            else
                C{FitNo,ExpNo} = num2str(C_values(FitNo,ExpNo));
                C_index(FitNo,ExpNo) = 0; 
            end
            if Tau_holds(ExpNo) == 0 % Lifetime
                if FitNo == 1
                    var = var + 1;
                    Tau_index(ExpNo) = var;
                end
                p_i         = num2str(Tau_index(ExpNo));
                p_start(Tau_index(ExpNo))= Tau_values(ExpNo);
                Tau{1,ExpNo}      = ['p(' p_i ')'];
                p_text{Tau_index(ExpNo)} = ['Tau' num2str(ExpNo)'];
                LB(Tau_index(ExpNo)) = 0; % Lifetimes can't be negative!
                UB(var) = +Inf;
            else
                Tau{1,ExpNo} = num2str(Tau_values(ExpNo));
            end
            if Beta_holds(ExpNo) == 0 % Beta
                if FitNo == 1
                    var = var + 1;
                    Beta_index(ExpNo) = var;
                end
                p_i         = num2str(Beta_index(ExpNo));
                p_start(Beta_index(ExpNo))= Beta_values(ExpNo);
                Beta{1,ExpNo}      = ['p(' p_i ')'];
                p_text{Beta_index(ExpNo)} = ['Beta' num2str(ExpNo)'];
                LB(Beta_index(ExpNo)) = 0; % 0 < Beta < 1
                UB(Beta_index(ExpNo)) = 1;
            else
                Beta{1,ExpNo} = num2str(Beta_values(ExpNo));
            end
            % Write the model
            if GaussIRF == 1
                expmodel_temp{FitNo,ExpNo}=['+0.5*',C{FitNo,ExpNo},'.*exp((',FWHM_f,'/(2*sqrt(2*log(2)))^2)/(2*(',Tau{1,ExpNo},')^2)).*exp(-(t-',t0_f,')/',Tau{1,ExpNo},').*(1+erf((t-',t0_f,'-(',FWHM_f,'/(2*sqrt(2*log(2)))^2)/',Tau{1,ExpNo},')/(sqrt(2)*',FWHM_f,'/(2*sqrt(2*log(2))))))'];
            elseif Stretched == 0
                expmodel_temp{FitNo,ExpNo}=['+' C{FitNo,ExpNo} '.*exp(-(t-' t0_f ')/' Tau{1,ExpNo} ').*heaviside(t-' t0_f ')'];
            elseif Stretched == 1
                expmodel_temp{FitNo,ExpNo}=['+' C{FitNo,ExpNo} '.*exp(-((t-' t0_f ')/' Tau{1,ExpNo} ').^' Beta{1,ExpNo} ').*heaviside(t-' t0_f ')'];
            end
        end      
    end
    for i=1:size(expmodel_temp,1)
        expmodel{i,1}={[expmodel_temp{i,:}]};
    end
end

% Now put together the pieces of the fit function
fitfunc_model=cell(nFits,1); % Preallocate
for FitNo=1:nFits
    if NumExp==0
        fitfunc_model{1}=infmodel{FitNo};
    else
        fitfunc_model{FitNo}=strcat(char(infmodel{FitNo}),char(expmodel{FitNo}));

    end
end
if fittype==1 % For local fits, give function handle
    fitfunc = eval(['@(p,t) ' fitfunc_model{1}]);
elseif fittype==2 % For global fits, give cell array
    fitfunc=cell(1,nFits);
    for m=1:nFits
        fitfunc{m} = eval(['@(p,t) ' fitfunc_model{m}]);
    end
end