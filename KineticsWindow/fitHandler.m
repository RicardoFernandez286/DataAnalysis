function [Dbest,res,CConc,Sfit,FitStruct] = fitHandler(doFit,PP_Data,KineticModel,ExtraPlots,Model_pFID,CutData,NormaliseSAS,Method)
%% Hardcoded settings
doParallel = true;
FillNaNs   = true;
%% Disable some warnings
warning('off','MATLAB:Axes:NegativeLimitsInLogAxis');
warning('off','MATLAB:Axes:NegativeDataInLogAxis');
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:illConditionedMatrix');
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:rankDeficientMatrix');
%% Extract IRF parameters and convert to fit parameters and boundaries
GauIRF      = KineticModel.IRF.GauIRF;  % Decide whether to do Gaussian IRF or not (otherwise Heaviside)
IRF_tbl     = cell2mat(KineticModel.IRF.Values(:,1:3));  % Numerical (IRF values)
IsFixIRF    = cell2mat(KineticModel.IRF.Values(:,4));    % Logical (fix or not)
IsT0        = [1 0];
if GauIRF == 1
    p0_IRF = IRF_tbl(~IsFixIRF,1)';
    UB_IRF = IRF_tbl(~IsFixIRF,2)';
    LB_IRF = IRF_tbl(~IsFixIRF,3)';
elseif GauIRF == 0
    p0_IRF = IRF_tbl(~IsFixIRF&IsT0,1)';
    UB_IRF = IRF_tbl(~IsFixIRF&IsT0,2)';
    LB_IRF = IRF_tbl(~IsFixIRF&IsT0,3)';
else
    error('Unknown parameter for Gaussian IRF configuration. Set it to 0 or 1 to disable/enable.')
end
%% Extract Kmat function and main parameters from kinetic model
Kmat_Fun    = KineticModel.Kmat_Fun;
numC        = KineticModel.numC;
numTaus     = KineticModel.numTaus;
C0          = KineticModel.C0;
modelType   = KineticModel.modelType;
%% Convert Rate table to fit parameters and boundaries
Tau0_tbl    = cell2mat(KineticModel.Tau0(:,1:3));  % Numerical (i.e. tau values)
IsFixTau    = cell2mat(KineticModel.Tau0(:,4));    % Logical (fix or not)
NfixTau     = sum(IsFixTau);
NfixIRF     = sum(IsFixIRF);
Nfix        = NfixTau + NfixIRF;
if Nfix == 0
    p0_Tau   = Tau0_tbl(:,1)';
    UB_Tau   = Tau0_tbl(:,2)';
    LB_Tau   = Tau0_tbl(:,3)';
    FixTAU   = [];
    FixIRF   = [];
else
    p0_Tau   = Tau0_tbl(~IsFixTau,1)';
    UB_Tau   = Tau0_tbl(~IsFixTau,2)';
    LB_Tau   = Tau0_tbl(~IsFixTau,3)';
    FixTAU   = Tau0_tbl(IsFixTau,1)';
    if GauIRF == 1
        FixIRF   = IRF_tbl(IsFixIRF,1)';
    else
        FixIRF   = IRF_tbl(IsFixIRF&IsT0,1)';
    end
end

% Infinite components have been taken care of (they are fixed).

%% Join all parameters into the global parameter vector
p0_all = [p0_IRF p0_Tau];
UB_all = [UB_IRF UB_Tau];
LB_all = [LB_IRF LB_Tau];
%% Read Pump--Probe dataset
Dexp   = PP_Data.rawsignal{1};
t      = PP_Data.delays;
WL     = PP_Data.cmprobe{1};
% Fill in missing NaNs
if sum(isnan(Dexp(:))) > 0 && FillNaNs == true
    Dexp = fillmissing(Dexp,'pchip');
    warning('Some NaN values were encountered in current dataset and were interpolated. If you do not want this behaviour, check the data and try again.')
end
%% Read Cut Settings and Cut Data
if CutData == 1
    tmin = PP_Data.plotranges{1}(1);
    tmax = PP_Data.plotranges{1}(2);
    Lmin = PP_Data.plotranges{1}(3);
    Lmax = PP_Data.plotranges{1}(4);
    selWL= WL>=Lmin & WL<=Lmax;
    selT = t>=tmin & t<=tmax;
    t = t(selT);
    WL= WL(selWL);
    Dexp=Dexp(selT,selWL);
end
%% Read Other Settings
linlog      = PP_Data.linlog;
timescale   = PP_Data.timescale;
probelabel  = PP_Data.probeunits;
Xunits      = PP_Data.Xunits;
%% Read Fit Options
p0eps = p0_all;
p0eps(p0_all==0)=1e-6; %#ok<*NASGU>
switch Method
    case {'Levenberg-Marquardt','Trust-Region-Reflective'}
        fitOptions = optimoptions('lsqnonlin',...
                        'MaxFunctionEvaluations',15000,...
                        'MaxIterations',350,...
                        'Algorithm',Method,...
                        'FiniteDifferenceType','central',...
                        'OptimalityTolerance',1e-15,...
                        'FunctionTolerance',1e-15,...
                        'StepTolerance',1e-15,...
                        'UseParallel',false,...
                        'ScaleProblem','Jacobian',...
                        'TypicalX',p0eps,...
                        'Display','off',...
                        'UseParallel',false,...
                        'PlotFcn',{'optimplotresnorm_log','optimplotx_txt'});
    case 'Particle Swarm'
        PSOoptions = optimoptions('particleswarm',...
                        'FunctionTolerance',1e-6,...
                        'MaxIterations',350,...
                        'Display','off',...
                        'PlotFcn',{'pswplotbestf','optimplot_bestx'},...
                        'UseVectorized',false,...
                        'UseParallel',doParallel);
        if sum(isinf(UB_all - LB_all))==0
            PSOoptions.InitialSwarmSpan = 2*(UB_all - LB_all);
        end
    case 'Genetic Algorithm'
        GAoptions = optimoptions('ga',...
            'Display','off',...
            'MaxGenerations',100,...
            'PlotFcn',{'gaplotbestf','gaplotbestindiv_txt'}, ...
            'UseParallel',doParallel, ...
            'UseVectorized',false);
    case 'Constrained Fit (fmincon)'
        FMCopt = optimoptions('fmincon', ...
            'ConstraintTolerance',1e-3, ...
            'FiniteDifferenceType','central', ...
            'MaxFunctionEvaluations',15000,...
            'MaxIterations',350,...
            'OptimalityTolerance',1e-15,...
            'StepTolerance',1e-15,...
            'Display','off',...
            'Algorithm','sqp',...
            'FunctionTolerance',1e-15,...
            'PlotFcn',{'optimplotfval','optimplotx_txt'});
           
            fprintf('\nNon-linear constraints require additional input:\n')
            t1 = input('  Enter Tau1 (TCSPC) in ns: ');
            t2 = input('  Enter Tau2 (TCSPC) in ns: ');
            a1 = input('  Enter A1   (TCSPC) in %:  ');
            a2 = 100-a1;
            TCSPC = [t1 t2 a1 a2];
end
%% Do fit or calculate
% Determine whether a fit is to be done or just plot the initial guess
if doFit == 1
    tic;
    Nvars = length(p0_all);
    switch Method
        case {'Levenberg-Marquardt','Trust-Region-Reflective'}
            [pFit,~,~,exitflag,output,~,J] = lsqnonlin(@(p) FitFunc(p,t,Kmat_Fun,C0,Dexp,GauIRF,FixTAU,IsFixTau,FixIRF,IsFixIRF,Model_pFID),p0_all,LB_all,UB_all,fitOptions);
        case 'Particle Swarm'
            fitfun = @(p) norm(FitFunc(p,t,Kmat_Fun,C0,Dexp,GauIRF,FixTAU,IsFixTau,FixIRF,IsFixIRF,Model_pFID)).^2;
           
            [pFit,~,exitflag,output] = particleswarm(fitfun,Nvars,LB_all,UB_all,PSOoptions);
        case 'Genetic Algorithm'
            fitfun = @(p) norm(FitFunc(p,t,Kmat_Fun,C0,Dexp,GauIRF,FixTAU,IsFixTau,FixIRF,IsFixIRF,Model_pFID)).^2;
           
            [pFit,~,exitflag,output] = ga(fitfun,Nvars,[],[],[],[],LB_all,UB_all,[],GAoptions);
        case 'Constrained Fit (fmincon)'
            fitfun = @(p) norm(FitFunc(p,t,Kmat_Fun,C0,Dexp,GauIRF,FixTAU,IsFixTau,FixIRF,IsFixIRF,Model_pFID)).^2;
            nlc = @(p) nlcon(p,TCSPC,FixTAU,IsFixTau,IsFixIRF);
           
            [pFit,~,exitflag,output,~,grad,hessian] = fmincon(fitfun,p0_all,[],[],[],[],LB_all,UB_all,nlc,FMCopt);
    end
   
    tfit=toc;
    [~,Dbest,CConc,Sfit] = FitFunc(pFit,t,Kmat_Fun,C0,Dexp,GauIRF,FixTAU,IsFixTau,FixIRF,IsFixIRF,Model_pFID);
    res     = Dexp-Dbest;
    chi2    = norm(res).^2./(numel(res)-length(pFit)-Nfix-1);
    if exitflag == -1
        pErr = NaN.*zeros(Nvars,1);
    else
        switch Method
            case {'Levenberg-Marquardt','Trust-Region-Reflective'}
                % First, calculate covariance matrix from the Jacobian
                covB = chi2.*inv(J'*J);
                % This calculates the 95 % confidence intervals
                ci   = nlparci(pFit,res,'covar',covB);
                % Take errors as the difference between the parameter and the lower ci
                pErr = diff(ci,1,2)./2;
            case 'Constrained Fit (fmincon)'
                % First, calculate covariance matrix from the Jacobian
                covB = chi2.*inv(hessian);
                % This calculates the 95 % confidence intervals
                ci   = nlparci(pFit,res,'covar',covB);
                % Take errors as the difference between the parameter and the lower ci
                pErr = diff(ci,1,2)./2;
            otherwise
                pErr = NaN.*zeros(Nvars,1);
        end
    end
else
    pFit    = p0_all;
    [~,Dbest,CConc,Sfit] = FitFunc(pFit,t,Kmat_Fun,C0,Dexp,GauIRF,FixTAU,IsFixTau,FixIRF,IsFixIRF,Model_pFID);
    res     = Dexp-Dbest;
    chi2    = norm(res).^2./(numel(res)-length(pFit)-Nfix-1);
    pErr    = NaN(size(pFit));
end
% Store in output structure (FitStruct) for further reference
FitStruct.pFit      = pFit;
FitStruct.pErr      = pErr;
FitStruct.GauIRF    = GauIRF;
FitStruct.FixTAU    = FixTAU;
FitStruct.FixIRF    = FixIRF;
FitStruct.IsFixTau  = IsFixTau;
FitStruct.IsFixIRF  = IsFixIRF;
FitStruct.chi2      = chi2;
FitStruct.FitDone   = doFit;
%% If the data was previously cut, fill with NaNs
if CutData == 1
    Dbest(selT,selWL)   = Dbest;
    Dexp(selT,selWL)    = Dexp;
    res(selT,selWL)     = res;
    CConc(selT,:)       = CConc;
    Sfit(:,selWL)       = Sfit;
    Dbest(~selT,:)      = NaN;
    Dbest(:,~selWL)     = NaN;
    Dexp(~selT,:)       = NaN;
    Dexp(:,~selWL)      = NaN;
    res(~selT,:)        = NaN;
    res(:,~selWL)       = NaN;
    CConc(~selT,:)      = NaN;
    Sfit(:,~selWL)      = NaN;
    t(selT)     = t;
    WL(selWL)   = WL;
    t(~selT)    = NaN;
    WL(~selWL)  = NaN;
end
%% Display fit results in command window
fprintf('\n**********************************');
fprintf('\n=======================\n\t Fit Results\n(%s)\n=======================',datetime)
fprintf('\nDataset: "%s"\n',PP_Data.datafilename)
switch modelType
    case 'Parallel'
        fprintf('\nFit Model: Parallel, %i components',numC);
    case 'Sequential'
        fprintf('\nFit Model: Sequential, %i components',numC);
    case 'Target'
        fprintf('\nFit Model: Target\n\t%s',KineticModel.TargetModel);
    case 'ODE-defined'
        fprintf('NOT IMPLEMENTED!');
end
fprintf('\n');
if doFit == 1
    fprintf('\nFit done in %.2f s\n',tfit);
    fprintf('Method: %s\n',Method);
    switch Method
        case {'Levenberg-Marquardt','Trust-Region-Reflective'}
            switch exitflag
                case 1
                    exit_str = 'Function converged to a solution x.';
                case 2
                    exit_str = 'Change in parameters is less than the specified tolerance, or Jacobian at current parameters is undefined.';
                case 3
                    exit_str = 'Change in the residual is less than the specified tolerance.';
                case 4
                    exit_str = 'Relative magnitude of search direction is smaller than the step tolerance.';
                case 0
                    exit_str = 'Number of iterations or function evaluations exceeded';
                case -1
                    exit_str = 'Optimization terminated by output or plot function.';
                case -2
                    exit_str = 'Problem is infeasible: the bounds lb and ub are inconsistent.';
            end
            fprintf('\t %i Iterations, %i Func. Evals.\n',output.iterations,output.funcCount);
        case 'Particle Swarm'
            switch exitflag
                case 1
                    exit_str = 'May have converged!, relative change in SSR over the last iterations is < FuncTol.';
                case 0
                    exit_str = 'Number of iterations exceeded.';
                case -1
                    exit_str = 'Optimization terminated by output or plot function.';
                case -2
                    exit_str = 'Problem is infeasible: the bounds LB and UB are inconsistent.';
                case -3
                    exit_str = 'The best objective function value is below the objective limit.';
                case -4
                    exit_str = 'Best objective function value did not change within allowed run time.';
                case -5
                    exit_str = 'Run time exceeded maximum allowed run time.';
            end
            fprintf('\t %i Iterations, %i Func. Evals.\n',output.iterations,output.funccount);
        case 'Genetic Algorithm'
            switch exitflag
                case 1
                    exit_str = 'Average cumulative change in SSR < FuncTol, and params. within bounds ± ConstTol.';
                case 3
                    exit_str = 'Change in SSR < FuncTol, and params. within bounds ± ConstTol.';
                case 4
                    exit_str = 'Magnitude of step smaller than eps and params. within bounds ± ConstTol.';
                case 5
                    exit_str = 'Maximum fitness reached.';
                case 0
                    exit_str = 'Maximum number of generations exceeded.';
                case -1
                    exit_str = 'Optimization terminated by output or plot function.';
                case -2
                    exit_str = 'No feasible point found.';
                case -4
                    exit_str = 'Stall time limit exceeded.';
                case -5
                    exit_str = 'Time limit exceeded.';
            end
        case 'Constrained Fit (fmincon)'
            switch exitflag
                case 1
                    exit_str = 'Gradient norm less than tolerance, params. within bounds ± ConstTol.';
                case 0
                    exit_str = 'Number of iterations exceeded.';
                case -1
                    exit_str = 'Optimization terminated by output or plot function.';
                case -2
                    exit_str = 'No feasible point found.';
                case 2
                    exit_str = 'Change in parameters less than StepTol, params. within bounds ± ConstTol.';
            end
            [~,nlc_val] = nlc(pFit);
            fprintf('\tNon-linear constraints: %.3g | %.3g | %.3g |\n',nlc_val);
            fprintf('\t%i Iterations, %i Func. Evals.\n',output.iterations,output.funcCount);
    end
    fprintf('\t[%s]\n\n',exit_str)
else
    fprintf('Preview only. No fit done.\n\n');
end
if IsFixIRF(1) == 1
    t0 = FixIRF(1);
    fprintf('t0   = %.3f %s**\n',t0,timescale)
else
    t0      = pFit(1);
    t0err   = pErr(1);
    fprintf('t0   = %.3f ± %.3f %s\n',t0,t0err,timescale)
end
if GauIRF == 1
    if IsFixIRF(2) == 1
        FWHM = FixIRF(end);
        fprintf('FWHM = %.3f %s**\n[Gaussian IRF assumed]\n\n',FWHM,timescale)
    else
        FWHM    = pFit(2-IsFixIRF(1));
        FWHMerr = pErr(2-IsFixIRF(1));
        fprintf('FWHM = %.3f ± %.3f %s\n[Gaussian IRF assumed]\n\n',FWHM,FWHMerr,timescale)
    end
else
    fprintf('[Delta-shaped IRF assumed]\n\n')
end
beginTau            = (1+length(IsFixIRF)-sum(IsFixIRF));
FitTaus             = pFit(beginTau:end);
ErrTaus             = pErr(beginTau:end);
TAU_fit             = zeros(size(Tau0_tbl(:,1)));
TAU_err             = zeros(size(Tau0_tbl(:,1)));
TAU_fit(IsFixTau)   = FixTAU;
TAU_fit(~IsFixTau)  = FitTaus;
TAU_err(IsFixTau)   = NaN;
TAU_err(~IsFixTau)  = ErrTaus;
Ntaus               = length(TAU_fit);
for i=1:Ntaus
    if IsFixTau(i)
        fprintf('tau%i = %.3f %s**\n',i,TAU_fit(i),timescale)
    else
        fprintf('tau%i = %.3f ± %.3f %s\n',i,TAU_fit(i),TAU_err(i),timescale)
    end
end
    fprintf('\n# of fixed IRF parameters: %i',sum(IsFixIRF));
    fprintf('\n# of fixed Time Constants: %i',sum(IsFixTau));
fprintf('\n\nchi^2 = %.4g\n',chi2);
fprintf('**********************************\n');
%% Plot results
%%% Plot the concentration profiles and SAS
% fhA     = figure(3);
% clf(fhA);
% fhA.Name= 'SAS and Temporal Evolution (Concentration Profiles)';
% fhA.Position(3:4) = [700 1060];
% movegui(fhA,'onscreen')
%
% axA = axes('parent',fhA);
FontSize = 16;
nC      = size(CConc,2);
if nC == 2
    cmapR= [0 0 0.85; 0.85 0 0];
else
    cmapR= brighten(turbo(nC+1),-0.5);
end
%%%%% Concentration profiles
% ax = subplot(2,1,1,axA);
fhA1     = figure(3);
clf(fhA1);
fhA1.Name= 'Temporal Evolution (Concentration Profiles)';
fhA1.Position(3:4) = [700 410];
movegui(fhA1,'onscreen')
ax = axes('parent',fhA1);
hold(ax,'on')
for i=1:nC
    plot(ax,t,CConc(:,i),'LineWidth',2,'Color',cmapR(i,:),'DisplayName',char(64+i));
end
hold(ax,'off')
title(ax,'Concentration Profiles');
xline(ax,t0);
yline(ax,0);
axis(ax,'tight');
ax.Box = 'on';
ax.FontSize = FontSize;
ax.XScale = linlog;
xlabel(ax,['Delay (' timescale ')'],'FontWeight','bold')
ylabel(ax,'Relative Conc.','FontWeight','bold')
if GauIRF
    if IsFixIRF(2) == 1
        FWHM  = FixIRF(end);
    else
        FWHM  = pFit(2);
    end
    xlim(ax,[round(10*FWHM)/10 max(t)]);
else
    xlim(ax,[min(t(t>0)) max(t)]);
end
ax.TickLength   = [0.02 0.02];
% ax.TickDir      = 'both';
%%%%%% SAS
% ax = subplot(2,1,2);
fhA2     = figure(4);
clf(fhA2);
fhA2.Position(3:4) = [700 410];
fhA2.Position(2) = fhA1.Position(2) - 500;
movegui(fhA2,'onscreen')
ax = axes('parent',fhA2);
switch modelType
    case 'Parallel'
        SpectraTitle = 'Decay-Associated Spectra';
    otherwise
        SpectraTitle = 'Species-Associated Spectra';
end
fhA2.Name= SpectraTitle;
hold(ax,'on')
% Automatically detect discontinuities in the probe axis (e.g. masked pump scatter)
% and plot them accordingly
dProbe          = diff(WL) - mean(diff(WL));
jump_ID         = dProbe >= 5*mean(diff(WL));
Sfit(:,jump_ID) = NaN;
for i=1:nC
    switch modelType
        case 'Parallel'
            % Convert Tau to a readable format (can handle tau < 1 s)
            TAU_plot = TAU_fit(i);
            powS     = getTimescale(timescale);
            newExp   = 3.*floor(floor(log10(TAU_plot))'./3);
            powS     = powS + newExp;
            TAU_plot = TAU_plot.*10.^(-newExp);
            ERR_plot = TAU_err(i).*10.^(-newExp);
            tS_plot  = setTimescale(powS);
           
            if isinf(TAU_fit(i))
                tS_plot = '';
            end
            if ~isinf(TAU_fit(i))
                name = ['\tau_{' char(64+i) '} = ' num2str(TAU_plot,'%4.3g')];
            else
                name = ['\tau_{' char(64+i) '} = \infty'];
            end
            if IsFixTau(i) && isinf(TAU_fit(i))
                name = [name '*']; %#ok<*AGROW>
            elseif IsFixTau(i) && ~isinf(TAU_fit(i))
                name = [name ' ' tS_plot '*'];
            else
                name = [name '\pm' num2str(round(ERR_plot,1,'significant'),'%2.2g') ' ' tS_plot];
            end
        case 'Sequential'
            % Convert Tau to a readable format (can handle tau < 1 s)
            TAU_plot = TAU_fit(i);
            powS     = getTimescale(timescale);
            newExp   = 3.*floor(floor(log10(TAU_plot))'./3);
            powS     = powS + newExp;
            TAU_plot = TAU_plot.*10.^(-newExp);
            ERR_plot = TAU_err(i).*10.^(-newExp);
            tS_plot  = setTimescale(powS);
            if ~isinf(TAU_fit(i))
                name = ['\tau_{' num2str(i) '} = ' num2str(TAU_plot,'%4.3g')];
            else
                name = ['\tau_{' num2str(i) '} = \infty'];
            end
            if IsFixTau(i) && isinf(TAU_fit(i))
                name = [name '*'];
            elseif IsFixTau(i) && ~isinf(TAU_fit(i))
                name = [name ' ' tS_plot '*'];
            else
                name = [name '\pm' num2str(round(ERR_plot,1,'significant'),'%2.2g') ' ' tS_plot];
            end
        otherwise
            name = char(64+i);
    end
    switch NormaliseSAS
        case 0
            plot(WL,Sfit(i,:),'Color',cmapR(i,:),'LineWidth',2,'DisplayName',name);
            Yname = 'Amplitude (mOD)';
        case 1
            plot(WL,Sfit(i,:)./max(abs(Sfit(i,:))),'Color',cmapR(i,:),'LineWidth',2,'DisplayName',name);
            Yname = 'Norm. SAS Amplitude';
    end
   
end
hold(ax,'off')
legend('location','best');
yline(ax,0,'handlevisibility','off');
title(ax,SpectraTitle);
xlim(ax,'tight');
ax.Box      = 'on';
ax.FontSize = FontSize;
xlabel(ax,[probelabel ' (' Xunits ')'],'FontWeight','bold')
ylabel(ax,Yname,'FontWeight','bold')
ax.TickLength   = [0.02 0.02];
% ax.TickDir      = 'both';
drawnow;
%%% Plot SVD of residuals
% [U,s,W] = svd(res(selT,selWL),'vector');
% plotSVD = 4;
% clf(figure(4));
% plot(t(selT),U(:,plotSVD)'); ax=gca; ax.XScale = 'log'; yline(0);
% plot(WL(selWL),W(:,plotSVD)'); yline(0);
if ExtraPlots == 0
    return
end
%% Contour plots of the dataset, fit and residuals
fhB     = figure(5);
clf(fhB);
fhB.Name= 'Data, Fit and Residuals (Contour Plots)';
fhB.Position(3:4) = [1430 420]; % Resize the figure to a wide shape
movegui(fhB,'onscreen')
figure(fhB);
tiledlayout(fhB,1,3,'padding','compact');
tplot   = PP_Data.delays-t0;    % Shift time axis such that tplot = 0 at t0
WLplot  = PP_Data.cmprobe{1};
zoomF   = 100/100;      % A zoom factor for the Z scale of the data
zoomR   = 10/100;       % A zoom factor for the Z scale of the residuals
NCtrs   = 40;           % Number of contours to plot
% Get minimum and maximum absorbance of the dataset
% maxZ    = max(Dexp(:));
% minZ    = min(Dexp(:));
maxZ    = max(abs(Dexp(:)));
minZ    = -maxZ;
% Get minimum and maximum absorbance of the residuals
maxZres = max(abs([maxZ minZ]));
minZres = -max(abs([maxZ minZ]));
% Define contours to plot
ctrs    = linspace(zoomF*minZ,zoomF*maxZ,NCtrs+1);
ctrs_r  = linspace(zoomR*minZres,zoomR*maxZres,NCtrs+1);
% Plot the data
ax1 = nexttile;
Dexp(:,jump_ID) = NaN;
contourf(ax1,WLplot,tplot,Dexp,ctrs,'EdgeColor','flat'); cb1 = colorbar;
colormap(ax1,darkb2r(zoomF*minZ,zoomF*maxZ,NCtrs,2));
caxis(ax1,zoomF*[minZ maxZ]);
ax1.YScale = linlog;
% Plot the fit
ax2 = nexttile;
Dbest(:,jump_ID) = NaN;
contourf(ax2,WLplot,tplot,Dbest,ctrs,'EdgeColor','flat'); cb2 = colorbar;
colormap(ax2,darkb2r(zoomF*minZ,zoomF*maxZ,NCtrs,2));
caxis(ax2,zoomF*[minZ maxZ]);
title(ax2,'Fit','FontSize',18);
ax2.YScale = linlog;
% Plot the residuals
ax3 = nexttile;
res(:,jump_ID) = NaN;
contourf(ax3,WLplot,tplot,res,ctrs_r,'EdgeColor','flat'); cb3 = colorbar;
colormap(ax3,darkb2r(zoomR*minZres,zoomR*maxZres,NCtrs,2));
caxis(ax3,zoomR*[minZres maxZres]);
ax3.YScale = linlog;
% Link the axes of the three plots
linkaxes([ax1,ax2,ax3],'xy');
% Add lines at time zero
yline(ax1,0,'-','Color',[0.5 0.5 0.5],'LineWidth',1.5);
yline(ax2,0,'-','Color',[0.5 0.5 0.5],'LineWidth',1.5);
yline(ax3,0,'-','Color',[0.5 0.5 0.5],'LineWidth',1.5);
% Add titles and axes labels
title(ax1,'(A) Data','FontAngle','italic');
title(ax2,'(B) Fit','FontAngle','italic');
title(ax3,'(C) Residuals','FontAngle','italic');
xlabel(ax1,[probelabel ' (' Xunits ')'],'FontWeight','bold')
xlabel(ax2,[probelabel ' (' Xunits ')'],'FontWeight','bold')
xlabel(ax3,[probelabel ' (' Xunits ')'],'FontWeight','bold')
ylabel(ax1,['Time (' timescale ')'],'FontWeight','bold')
ylabel(ax2,['Time (' timescale ')'],'FontWeight','bold')
ylabel(ax3,['Time (' timescale ')'],'FontWeight','bold')
% Increase Font Size
ax1.FontSize = 14;
ax2.FontSize = 14;
ax3.FontSize = 14;
% Show ticks inside and outside and make them longer
ax1.TickDir = 'both';
ax2.TickDir = 'both';
ax3.TickDir = 'both';
ax1.TickLength = [0.02 0.02];
ax2.TickLength = [0.02 0.02];
ax3.TickLength = [0.02 0.02];
% Add title to the colorbars
title(cb1,{'\DeltaAbs (mOD)'})
title(cb2,{'\DeltaAbs (mOD)'})
title(cb3,{'\DeltaAbs (mOD)'})
% Change X+Y limits to the fitted ranges
xlim(ax1,[Lmin Lmax]);
if tmin < 0
    tminPlot = min(tplot(tplot>0));
else
    tminPlot = tmin;
end
ylim(ax1,[tminPlot tmax]);
drawnow;
%%%%%% END OF FILE
end
%% To-do:
% Model_pFID
% Ask about initial concentrations
%%%%%%%%%%%%%%%%% EOF
%% This auxilliary function calculates the fit dataset
function [res,Dfit,CConc,Sfit] = FitFunc(p,t,Kmodel,C0,Dexp,GauIRF,FixTAU,IsFixTau,FixIRF,IsFixIRF,Model_pFID)
    beginTau            = (1+length(IsFixIRF)-sum(IsFixIRF));
    FitTaus             = p(beginTau:end);
    AllRates(IsFixTau)  = 1./FixTAU;
    AllRates(~IsFixTau) = 1./FitTaus;
    % Read parameters off the parameter vector (p)
    if IsFixIRF(1) == 1
        t0  = FixIRF(1);
    else
        t0  = p(1);
    end
    % Build the K matrix. Note that in the fit we are using k=1/tau
    Kmat    = Kmodel(AllRates);
    Conc_fn = @(t) kineticsKmat_simu(t,Kmat,C0,0,[]);
    if GauIRF == 1
        if IsFixIRF(2) == 1
            FWHM  = FixIRF(end);
        else
            FWHM  = p(2);
        end
        CConc   = IRF_Kmat(t,Kmat,C0,FWHM,t0);
%         IRFfnc  = @(T,T0,W) (2./W)*sqrt(log(2)/pi)*exp(-4*log(2).*((T-T0)./W).^2)';
%         CConc   = IRFconvol(IRFfnc,t,t0,FWHM,Conc_fn);
    else
        % Heaviside step function, usual definition with H(0) = 1/2
        Hfunc   = @(t) (1+sign(t))/2;
        % Multiply concentration profiles by Heaviside step function
        CConc = Conc_fn(t).*Hfunc(t-t0);
    end
    % Calculate Sfit and Dfit matrices (linear least squares)
    Sfit    = CConc\Dexp;
    Dfit    = CConc*Sfit;
   
    % Calculate Residuals
    res     = Dexp-Dfit;
end
%% This auxiliary function calculates the nonlinear constrains for a specific model
function [c,ceq] = nlcon1(p,TCSPC,FixTAU,IsFixTau,IsFixIRF)
    beginTau  = (1+length(IsFixIRF)-sum(IsFixIRF));
    FitTaus             = p(beginTau:end);
    AllTaus(IsFixTau)   = FixTAU;
    AllTaus(~IsFixTau)  = FitTaus;
    k         = 1./AllTaus;
   
    c = [];
    lam1 = 0.5.*(-sum(k(2:5)) - sqrt(sum(k(2:5)).^2 - 4.*(k(3).*k(4)+k(2).*k(5)+k(3).*k(5)))) + 1./TCSPC(1);
    lam2 = 0.5.*(-sum(k(2:5)) + sqrt(sum(k(2:5)).^2 - 4.*(k(3).*k(4)+k(2).*k(5)+k(3).*k(5)))) + 1./TCSPC(2);
   
    ceq = [lam1,lam2];
end
%% Auxiliary function, constraint 2
function [c,ceq] = nlcon(p,TCSPC,FixTAU,IsFixTau,IsFixIRF)
    beginTau  = (1+length(IsFixIRF)-sum(IsFixIRF));
    FitTaus             = p(beginTau:end);
    AllTaus(IsFixTau)   = FixTAU;
    AllTaus(~IsFixTau)  = FitTaus;
    k         = 1./AllTaus;
   
    c = [];
    lam1 = 0.5.*(-sum(k(2:5)) - sqrt(sum(k(2:5)).^2 - 4.*(k(3).*k(4)+k(2).*k(5)+k(3).*k(5)))) + 1./TCSPC(1);
    lam2 = 0.5.*(-sum(k(2:5)) + sqrt(sum(k(2:5)).^2 - 4.*(k(3).*k(4)+k(2).*k(5)+k(3).*k(5)))) + 1./TCSPC(2);
    theRatio = (k(4)+k(5)-(1./TCSPC(2)))./((1./TCSPC(1))-k(4)-k(5)) - TCSPC(4)./TCSPC(3);
   
    ceq = [lam1,lam2,theRatio];
end