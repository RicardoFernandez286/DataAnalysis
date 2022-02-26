function [Dbest,res,CConc,Sfit] = fitHandler(doFit,PP_Data,KineticModel,ExtraPlots,Model_pFID)

%% Disable some warnings
warning('off','MATLAB:Axes:NegativeLimitsInLogAxis');
warning('off','MATLAB:Axes:NegativeDataInLogAxis');
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:illConditionedMatrix');
warning('off','MATLAB:nearlySingularMatrix');

%% Extract IRF parameters and convert to fit parameters and boundaries
GauIRF  = KineticModel.IRF.GauIRF;  % Decide whether to do Gaussian IRF or not (otherwise Heaviside)
t0_param= KineticModel.IRF.t0;      % Val UB LB
w_param = KineticModel.IRF.FWHM;    % Val UB LB

if GauIRF == 1
    p0_IRF = [t0_param(1) w_param(1)];
    UB_IRF = [t0_param(2) w_param(2)];
    LB_IRF = [t0_param(3) w_param(3)];
else
    p0_IRF = t0_param(1);
    UB_IRF = t0_param(2);
    LB_IRF = t0_param(3);
end

%% Extract Kmat function and main parameters from kinetic model
Kmat_Fun = KineticModel.Kmat_Fun;
% numC     = KineticModel.numC;
% numTaus  = KineticModel.numTaus;
C0       = KineticModel.C0;
modelType= KineticModel.modelType;
%% Convert Rate table to fit parameters and boundaries
Tau0_tbl = cell2mat(KineticModel.Tau0(:,1:3));  % Numerical (i.e. tau values)
IsFixTau = cell2mat(KineticModel.Tau0(:,4));    % Logical (fix or not)
Nfix     = sum(IsFixTau);

if Nfix == 0
    p0_Tau   = Tau0_tbl(:,1)';
    UB_Tau   = Tau0_tbl(:,2)';
    LB_Tau   = Tau0_tbl(:,3)';
    FixP     = [];
else
    p0_Tau   = Tau0_tbl(~IsFixTau,1)';
    UB_Tau   = Tau0_tbl(~IsFixTau,2)';
    LB_Tau   = Tau0_tbl(~IsFixTau,3)';
    FixP     = Tau0_tbl(IsFixTau,1)';
end

% If we have any infinite components, their lower bounds 
%% Join all parameters into the global parameter vector
p0_all = [p0_IRF p0_Tau];
UB_all = [UB_IRF UB_Tau];
LB_all = [LB_IRF LB_Tau];

%% Read Pump--Probe dataset
Dexp   = PP_Data.rawsignal{1};
t      = PP_Data.delays;
WL     = PP_Data.cmprobe{1};
% WL = 1:length(WL);
linlog      = PP_Data.linlog;
timescale   = PP_Data.timescale;
probelabel  = PP_Data.probeunits;
Xunits      = PP_Data.Xunits;

%% Now we are ready to set the fit options (TO BE READ IN THE FUTURE, CONTINUEHERE)
p0eps = p0_all;
p0eps(p0_all==0)=eps;

fitOptions = optimoptions('lsqnonlin',...
                'MaxFunctionEvaluations',15000,...
                'MaxIterations',350,...
                'Algorithm','levenberg-marquardt',...
                'FiniteDifferenceType','central',...
                'OptimalityTolerance',1e-15,...
                'FunctionTolerance',1e-15,...
                'StepTolerance',1e-15,...
                'UseParallel',false,...
                'ScaleProblem','Jacobian',...
                'Display','off',...
                'PlotFcn',{'optimplotresnorm_log','optimplotx'});
    
% Determine whether a fit is to be done or just plot the initial guess
if doFit == 1
    [pFit,resnorm,~,exitflag,output,~,J] = lsqnonlin(@(p) FitFunc(p,t,Kmat_Fun,C0,Dexp,GauIRF,FixP,IsFixTau),p0_all,LB_all,UB_all,fitOptions);
    [~,Dbest,CConc,Sfit] = FitFunc(pFit,t,Kmat_Fun,C0,Dexp,GauIRF,FixP,IsFixTau);
    res     = Dexp-Dbest;
    chi2    = resnorm./(numel(res)-length(pFit)-Nfix-1);
else
    pFit    = p0_all;
    [~,Dbest,CConc,Sfit] = FitFunc(pFit,t,Kmat_Fun,C0,Dexp,GauIRF,FixP,IsFixTau);
    res     = Dexp-Dbest;
    chi2    = NaN;
end

%% Display fit results in command window
fprintf('\n=============\nFit Results:\n=============\n')
fprintf('t0   = %.3f \n',pFit(1))
fprintf('FWHM = %.3f \n\n',pFit(2))

for i=3:length(pFit)
    fprintf('k%i = %.3e ; tau = %.3f\n',i-2,1./pFit(i),pFit(i))
end
fprintf('\nchi^2 = %.3g\n',chi2);
fprintf('*****************\n');

%% Plot results
drawnow;
%%% Plot the concentration profiles and SAS
fhA     = figure(3);
clf(fhA);
fhA.Name= 'SAS and Temporal Evolution (Concentration Profiles)';
fhA.Position(3:4) = [700 1060];
movegui(fhA,'onscreen')

axA = axes('parent',fhA);

FontSize = 16;

nC      = size(CConc,2);
if nC == 2
    cmapR= [0 0 0.85; 0.85 0 0];
else
    cmapR= brighten(turbo(nC+1),-0.5);
end


% Concentration profiles
ax = subplot(2,1,1,axA);
hold(ax,'on')
for i=1:nC
    plot(ax,t,CConc(:,i),'LineWidth',2,'Color',cmapR(i,:),'DisplayName',char(64+i));
end
hold(ax,'off')
title(ax,'Concentration Profiles (a.u.)');
xline(ax,pFit(1));
yline(ax,0);
axis(ax,'tight');
ax.Box = 'on';
ax.FontSize = FontSize;
ax.XScale = linlog;
xlabel(ax,['Delay (' timescale ')'],'FontWeight','bold')
ylabel(ax,'Relative Conc.','FontWeight','bold')

if GauIRF
    xlim(ax,[round(10*pFit(2))/10 max(t)]);
else
    xlim(ax,[min(t(t>0)) max(t)]);
end
ax.TickLength   = [0.02 0.02];
% ax.TickDir      = 'both';

% SAS
ax = subplot(2,1,2);
switch modelType
    case 'Parallel'
        SpectraTitle = 'Decay-Associated Spectra';
    otherwise
        SpectraTitle = 'Species-Associated Spectra';
end

hold(ax,'on')

FitTaus = pFit((2+GauIRF):end);
TAU_fit(IsFixTau)  = FixP;
TAU_fit(~IsFixTau) = FitTaus;

for i=1:nC
    % Convert Tau to a readable format (can handle tau < 1 s)
    TAU_plot = TAU_fit(i);
    powS     = getTimescale(timescale);
    newExp   = 3.*floor(floor(log10(TAU_plot))'./3);
    powS     = powS + newExp;
    TAU_plot = TAU_plot.*10.^(-newExp);
    tS_plot  = setTimescale(powS);

    switch modelType
        case 'Parallel'
            if ~isinf(TAU_fit(i))
                name = ['\tau_{' char(64+i) '} = ' num2str(TAU_plot,'%4.3g') ' ' tS_plot];
            else
                name = ['\tau_{' char(64+i) '} = \infty'];
            end
            if IsFixTau(i)
                name = [name '*']; %#ok<*AGROW> 
            end
        case 'Sequential'
            if ~isinf(TAU_fit(i))
                name = ['\tau_{' num2str(i) '} = ' num2str(TAU_plot,'%4.3g') ' ' tS_plot];
            else
                name = ['\tau_{' num2str(i) '} = \infty'];
            end
            if IsFixTau(i)
                name = [name '*'];
            end
        otherwise
            name = char(64+i);
    end
   
    plot(WL,Sfit(i,:),'Color',cmapR(i,:),'LineWidth',2,'DisplayName',name);
end
hold(ax,'off')
legend('location','best');
yline(ax,0,'handlevisibility','off');
title(ax,SpectraTitle);
xlim(ax,'tight');
ax.Box = 'on';
ax.FontSize = FontSize;

xlabel(ax,[probelabel '(' Xunits ')'],'FontWeight','bold')
ylabel(ax,'Amplitude (mOD)','FontWeight','bold')

ax.TickLength   = [0.02 0.02];
% ax.TickDir      = 'both';

drawnow;

if ExtraPlots == 0
    return
end

%% Contour plots of the dataset, fit and residuals
fhB     = figure(4);
clf(fhB);
fhB.Name= 'Data, Fit and Residuals (Contour Plots)';
fhB.Position(3:4) = [1430 420]; % Resize the figure to a wide shape
movegui(fhB,'onscreen')

figure(fhB);
tiledlayout(fhB,1,3,'padding','compact');

tplot   = t-pFit(1);    % Shift time axis such that tplot = 0 at t0

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
contourf(ax1,WL,tplot,Dexp,ctrs,'EdgeColor','flat'); cb1 = colorbar;
colormap(ax1,darkb2r(zoomF*minZ,zoomF*maxZ,NCtrs,2));
caxis(ax1,zoomF*[minZ maxZ]);
ax1.YScale = linlog;

% Plot the fit
ax2 = nexttile;
contourf(ax2,WL,tplot,Dbest,ctrs,'EdgeColor','flat'); cb2 = colorbar;
colormap(ax2,darkb2r(zoomF*minZ,zoomF*maxZ,NCtrs,2));
caxis(ax2,zoomF*[minZ maxZ]);
title(ax2,'Fit','FontSize',18);
ax2.YScale = linlog;

% Plot the residuals
ax3 = nexttile;
contourf(ax3,WL,tplot,res,ctrs_r,'EdgeColor','flat'); cb3 = colorbar;
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

xlabel(ax1,[probelabel '(' Xunits ')'],'FontWeight','bold')
xlabel(ax2,[probelabel '(' Xunits ')'],'FontWeight','bold')
xlabel(ax3,[probelabel '(' Xunits ')'],'FontWeight','bold')

ylabel(ax1,'Time (s)','FontWeight','bold')
ylabel(ax2,'Time (s)','FontWeight','bold')
ylabel(ax3,'Time (s)','FontWeight','bold')

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
title(cb1,'Abs.','FontWeight','bold')
title(cb2,'Abs.','FontWeight','bold')
title(cb3,'Abs.','FontWeight','bold')

drawnow;

%%%%%% END OF FILE
end

%% To-do:
% Model_pFID
% Ask about initial concentrations


%%%%%%%%%%%%%%%%% EOF

%% This auxilliary function calculates the fit dataset
function [res,Dfit,CConc,Sfit] = FitFunc(p,t,Kmodel,C0,Dexp,GauIRF,FixP,IsFixTau)
%     numTaus = length(IsFixTau);
    FitTaus = p((2+GauIRF):end);
    AllTaus(IsFixTau)  = 1./FixP;
    AllTaus(~IsFixTau) = 1./FitTaus;

    % Read parameters off the parameter vector (p)
    t0      = p(1);
    % Build the K matrix. Note that in the fit we are using k=1/tau
    Kmat    = Kmodel(AllTaus);
    Conc_fn = @(t) kineticsKmat_simu(t,Kmat,C0,0,[]);

    if GauIRF == 1
        FWHM    = p(2);
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