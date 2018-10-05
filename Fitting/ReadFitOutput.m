function handles = ReadFitOutput(parameters,handles,fitP,fitP_SD)
% This function will build the fit function used to fit the data
% The input is a structure with the following fields:
%         Tau_values    Tau1 Tau2 Tau3 Tau4
%         Tau_holds     Tau1 Tau2 Tau3 Tau4
%         C_values      c1   c2   c3   c4   cInf
%         C_holds       c1   c2   c3   c4   cInf
%         Param_values  FWHM t0   Offset
%         Param_holds   FWHM t0   Offset
%         DoFit         0=Cancel, 1=Fit, 2=Preview
%         NumExp        Number of exponentials to fit
%         Switches      GaussIRF  StretchedExp (1 or 0)   
%
% The equation to fit is the result of an offset (c0) added to the 
% convolution of a sum of NumExp and an infinite component
%
%   Ricardo Fernández-Terán
%   v1.2 / 2018-10-05

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
%% Build the model
% Var is a number corresponding to the variable passed later to the fit
% model
var=0;

% p_start is a vector containing the starting values of the parameters
p_start=[];

% p_text is a cell array containing all parameters' names in order
out_text={};

% Check whether FWHM, t0 and Offset are held, then set them as a fixed
% number or as a variable to fit
i=1;
if GaussIRF ~= 0
    if Param_holds(1) == 0 % FWHM
        var         = var+1;
        FWHM_f      = num2str(fitP(var),2);
        FWHM_SD     = num2str(fitP_SD(var),2);
        out_text{i} = ['FWHM = ' FWHM_f ' ' char(177) ' ' FWHM_SD ' ' handles.timescale];
    else
        FWHM_f = [num2str(Param_values(1)) ' ' handles.timescale ' (H)'];
        out_text{i} = ['FWHM = ' FWHM_f];
    end
end

i=2;
if Param_holds(2) == 0 % t0
    var         = var+1;
    t0_f        = num2str(fitP(var),2);
    t0_SD       = num2str(fitP_SD(var),2);
    out_text{i} = [sprintf('t\x2080') ' = ' t0_f ' ' char(177) ' ' t0_SD ' ' handles.timescale];
else
    t0_f = [num2str(Param_values(2)) ' ' handles.timescale ' (H)'];
    out_text{i} = [sprintf('t\x2080') ' = ' t0_f];
end

i=3;
if C_holds(6) == 0 % Offset
    var           = var+1;
    Offset_f      = num2str(fitP(var),2);
    Offset_SD     = num2str(fitP_SD(var),2);
    out_text{i}   = ['Offset = ' Offset_f ' ' char(177) ' ' Offset_SD ' mOD'];
else
    Offset_f = [num2str(C_values(6)) ' mOD (H)'];
    out_text{i} = ['Offset = ' Offset_f];
end

i=4;
% Put in the infinite component
if C_holds(5) == 0 % Infinite component
    var         = var+1;
    cInf_f      = num2str(fitP(var),2);
    cInf_SD     = num2str(fitP_SD(var),2);
    out_text{i} = ['Infinite = ' cInf_f ' ' char(177) ' ' cInf_SD ' mOD'];
else
    cInf_f = [num2str(C_values(5)) ' mOD (H)'];
    out_text{i} = ['Infinite = ' cInf_f];
end
i=i+1;
out_text{i} = ' ';

% Build the exponentials
maxExp=4; % Set the maximum number of exponentials - for future upgrades
if NumExp ~= 0
    % Preallocate memory and initialise variables
    C = {}; C_SD = {};
    Beta = {}; Beta_SD = {};
    Tau = {}; Tau_SD = {}; TauAvg = {};
    % Check each exponential and write the model
    for k=1:NumExp
        i=i+1;
        if C_holds(k) == 0 % Amplitude
            var         = var+1;
            C{k}        = num2str(fitP(var),2);
            C_SD{k}     = num2str(fitP_SD(var),2);
            out_text{i+1} = ['C' sprintf(['\x' num2str(2080+k)]) ' = ' C{k} ' ' char(177) ' ' C_SD{k} ' mOD'];
        else
            C{k} = num2str(C_values(k));
            out_text{i+1} = ['C' sprintf(['\x' num2str(2080+k)]) ' = ' C{k} ' mOD (H)'];
        end
        i=i+1;
        if Tau_holds(k) == 0 % Lifetime
            var         = var+1;
            TauNumber(k)=fitP(var);
            if fitP(var) < 1000
                Tau{k}      = num2str(fitP(var),'%.3g');
                Tau_SD{k}   = num2str(fitP_SD(var),'%.3g');
                                    % TAU                    _k
                out_text{i-1} = [sprintf('\x3c4') sprintf(['\x' num2str(2080+k)]) ' = ' Tau{k} ' ' char(177) ' ' Tau_SD{k} ' ' handles.timescale];
            else
                Tau{k}      = num2str(fitP(var)/1000,2);
                Tau_SD{k}   = num2str(fitP_SD(var)/1000,2);
                switch handles.timescale
                    case 'ns'
                        longtimescale = [char(956) 's'];
                    case 'ps'
                        longtimescale = 'ns';
                    case 'fs'
                        longtimescale = 'ps';
                end
                out_text{i-1} = [sprintf('\x3c4') sprintf(['\x' num2str(2080+k)]) ' = ' Tau{k} ' ' char(177) ' ' Tau_SD{k} ' ' longtimescale];
            end
        else
            TauNumber(k)=Tau_values(k);
            if Tau_values(k) < 1000
                Tau{k} = num2str(Tau_values(k));
                out_text{i-1} = [sprintf('\x3c4') sprintf(['\x' num2str(2080+k)]) ' = ' Tau{k} ' ' handles.timescale ' (H)'];
            else
                Tau{k} = num2str(Tau_values(k)/1000);
                switch handles.timescale
                    case 'ns'
                        longtimescale = [char(956) 's'];
                    case 'ps'
                        longtimescale = 'ns';
                    case 'fs'
                        longtimescale = 'ps';
                end
                out_text{i-1} = [sprintf('\x3c4') sprintf(['\x' num2str(2080+k)]) ' = ' Tau{k} ' ' longtimescale ' (H)'];
            end
        end
        if Stretched == 1
            i=i+1;
            if Beta_holds(k) == 0 % Beta
                var             = var+1;
                Beta{k}         = num2str(fitP(var),2);
                Beta_SD{k}      = num2str(fitP_SD(var),2);
                BetaNumber(k)=fitP(var);
                out_text{i}   = [sprintf('\x3b2') sprintf(['\x' num2str(2080+k)]) ' = ' Beta{k} ' ' char(177) ' ' Beta_SD{k}];
            else
                Beta{k} = num2str(Beta_values(k));
                BetaNumber(k)=Beta_values(k);
                out_text{i} = [sprintf('\x3b2') sprintf(['\x' num2str(2080+k)]) ' = ' Beta{k} ' (H)'];
            end
            TauAvgNumber(k) = (TauNumber(k)/BetaNumber(k))*gamma(1/BetaNumber(k));
            if TauAvgNumber >= 1000
                TauAvg{k} = num2str(TauAvgNumber(k));
                switch handles.timescale   
                   case 'ns'
                       longtimescale = [char(956) 's'];
                   case 'ps'
                       longtimescale = 'ns';
                   case 'fs'
                       longtimescale = 'ps';
                end
            else
                longtimescale = handles.timescale;
                TauAvg{k} = num2str(TauAvgNumber(k),'%.3g');
            end
            out_text{i+1} = ['<' sprintf('\x3c4') sprintf(['\x' num2str(2080+k)]) '> = ' TauAvg{k} ' ' longtimescale];    
        end
        % Uncomment this for an extra line between each Tau component (cosmetics)
        i=i+1+Stretched;
        out_text{i} = ' ';
    end
end
handles.FitOutput_text = transpose(out_text);