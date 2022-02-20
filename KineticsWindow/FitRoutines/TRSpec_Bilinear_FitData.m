% Description:  This script fits time-dependent spectra based on coupled first-order reactions.
%               The K matrix model can be defined as a function, which takes as parameters a vector with the rate constants.
%               It can be easily adapted to more complex and arbitrary kinetic schemes represented by systems of ODEs.
%
% Tested and implemented in MATLAB R2021b
% v1.0

%% Hardcoded Settings and options
warning('off','MATLAB:rankDeficientMatrix') % Disable warnings about rank-defficient matrices.

doFit = 1;  % Decide whether to plot the initial guess (=0) or to do the fit (=1).
reFit = 0;  % Decide whether to use the previously fitted parameters as initial guess to the fit (=1) or not (=0).

%% Load Data
[fileName,filePath] = uigetfile({'*.csv; *.dat; *.pdat','ASCII Data Files'},'MultiSelect','off');
data    = readmatrix([filePath filesep fileName]);

if fileName == 0
    error('No file selected. Aborting.');
end

% Read time and wavelength vectors, and absorbance matrix (Dexp)
t       = data(1,2:end)';
WL      = data(2:end,1);
Dexp    = data(2:end,2:end)';
Dexp    = fillmissing(Dexp, 'linear');

%%% Cut after a certain time delay [uncomment if needed]
% tCut    = 600;
% Dexp    = Dexp(t<tCut,:);
% t       = t(t<tCut);

%% Initialise fit parameters and define model
% Initial concentrations vector
C0      = [1 0 0];
% Rate constants
k0      = [3 65]; % Example starting parameters for glucose [k1 k2]

%%% Define the model here, by editing the shape of the K matrix
%   The initial concentration vector (C0) and initial rate constants (k0) need to be updated accordingly
%   The model exemplified here is A-(k1)->B-(k2)->C-(k3)->D
model   = length(k0)+1; % Select model by number of rates
switch model
    case 4 % 4 components
        Kmodel = @(k) [-k(1)    0           0       0
                       +k(1)    -k(2)       0       0
                       0        +k(2)       -k(3)   0
                       0        0           +k(3)   0];
    case 3 % 3 components
        Kmodel = @(k) [-k(1)    0       0
                       +k(1)    -k(2)   0
                       0        +k(2)   0];
end

% Set the upper and lower bounds for the fitted rates
LB_k    = 1e-3*ones(1,length(k0));
UB_k    = 2000*ones(1,length(k0));

%%% Define the IRF parameters (t0 and FWHM, respectively)
IRF_0   = [2    1]; % Example starting parameters for glucose
LB_IRF  = [0    1];
UB_IRF  = [3    2];

% Join all parameters and bounds into the parameter vector
p0  = [IRF_0 k0];
LB  = [LB_IRF LB_k];
UB  = [UB_IRF UB_k];

%% Do the fit
p0(p0==0)=eps; % add a small number in case any of the parameters is identically zero
options = optimoptions('lsqnonlin',...
                    'MaxFunctionEvaluations',15000,...
                    'MaxIterations',200,...
                    'Algorithm','Levenberg-Marquardt',...
                    'OptimalityTolerance',1e-3,...
                    'TypicalX',p0,...
                    'FunctionTolerance',1e-3,...
                    'StepTolerance',1e-5,...
                    'UseParallel',false,...
                    'Display','final',...
                    'Display','off',...
                    'PlotFcn',{'optimplotresnorm','optimplotx'});

% Determine whether we are re-fitting a dataset
if reFit == 1 && exist('pFit','var') ~= 0
    p0 = pFit;
end

% Determine whether a fit is to be done or just plot the initial guess
if doFit == 1
    [pFit,resnorm,~,~,output,~,J] = lsqnonlin(@(p) FitFunc(p,t,Kmodel,C0,Dexp),p0,LB,UB,options);
    [~,Dbest,CConc,Sfit] = FitFunc(pFit,t,Kmodel,C0,Dexp);
    res     = Dexp-Dbest;
    chi2    = resnorm./(numel(res)-length(pFit)-1);
else
    pFit    = p0;
    chi2    = NaN;
    [~,Dbest,CConc,Sfit] = FitFunc(pFit,t,Kmodel,C0,Dexp);
    res     = Dexp-Dbest;
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
%%% Plot the concentration profiles and SAS
fhA     = figure(1);
clf(fhA);
fhA.Position(3:4) = [599 820];
axA = axes('parent',fhA);

nC      = size(CConc,2);
cmapR   = brighten(turbo(nC+1),-0.25);

% Concentration profiles
ax = subplot(2,1,1,axA);
hold(ax,'on')
for i=1:nC
    plot(ax,t-pFit(1),CConc(:,i),'LineWidth',2,'Color',cmapR(i,:),'DisplayName',['Species ' num2str(i)]);
end
hold(ax,'off')
title(ax,'Concentration Profiles (a.u.)');
xline(ax,0*pFit(1));
yline(ax,0);
axis(ax,'tight');
ax.Box = 'on';
ax.FontSize = 16;

% SAS
ax = subplot(2,1,2);
hold(ax,'on')
for i=1:nC
    plot(WL,Sfit(i,:),'Color',cmapR(i,:),'LineWidth',2,'DisplayName',['Species ' num2str(i)]);
end
hold(ax,'off')
legend;
yline(ax,0,'handlevisibility','off');
title(ax,'SAS');
xlim(ax,[400,750]);
ax.Box = 'on';
ax.FontSize = 16;

%%% Contour plots of the dataset, fit and residuals
fhB     = figure(2);
clf(fhB);
tiledlayout(fhB,1,3,'padding','compact');

tplot   = t-pFit(1);    % Shift time axis such that tplot = 0 at t0

zoomF   = 100/100;      % A zoom factor for the Z scale of the data
zoomR   = 10/100;       % A zoom factor for the Z scale of the residuals
NCtrs   = 40;           % Number of contours to plot

% Get minimum and maximum absorbance of the dataset
maxZ    = max(Dexp(:));
minZ    = min(Dexp(:));

% Get minimum and maximum absorbance of the residuals
maxZres = max(abs([maxZ minZ]));
minZres = -max(abs([maxZ minZ]));

% Define contours to plot
ctrs    = linspace(zoomF*minZ,zoomF*maxZ,NCtrs+1);
ctrs_r  = linspace(zoomR*minZres,zoomR*maxZres,NCtrs+1);

% Plot the data
ax1 = nexttile;
contourf(ax1,WL,tplot,Dexp,ctrs,'EdgeColor','flat'); cb1 = colorbar;
colormap(ax1,turbo(NCtrs)); 
caxis(ax1,zoomF*[minZ maxZ]);

% Plot the fit
ax2 = nexttile;
contourf(ax2,WL,tplot,Dbest,ctrs,'EdgeColor','flat'); cb2 = colorbar;
colormap(ax2,turbo(NCtrs));
caxis(ax2,zoomF*[minZ maxZ]);
title(ax2,'Fit','FontSize',18);

% Plot the residuals
ax3 = nexttile;
contourf(ax3,WL,tplot,res,ctrs_r,'EdgeColor','flat'); cb3 = colorbar;
colormap(ax3,darkb2r(zoomR*minZres,zoomR*maxZres,NCtrs,2));
caxis(ax3,zoomR*[minZres maxZres]);

% Link the axes of the three plots
linkaxes([ax1,ax2,ax3],'xy');

% Add lines at time zero
yline(ax1,0,'--','Color','w','LineWidth',2);
yline(ax2,0,'--','Color','w','LineWidth',2);
yline(ax3,0,'--','Color',[0.5 0.5 0.5],'LineWidth',2);

% Add titles and axes labels
title(ax1,'(A) Data','FontAngle','italic');
title(ax2,'(B) Fit','FontAngle','italic');
title(ax3,'(C) Residuals','FontAngle','italic');

xlabel(ax1,'Wavelength (nm)','FontWeight','bold')
xlabel(ax2,'Wavelength (nm)','FontWeight','bold')
xlabel(ax3,'Wavelength (nm)','FontWeight','bold')

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

% Change the Y scale (time) to logarithmic [uncomment if needed]
% ax1.YScale = 'log';
% ax2.YScale = 'log';
% ax3.YScale = 'log';

% Resize the figure to a wide shape
fhB.Position(3:4) = [1430 420];

%%%%%% END OF FILE
%%%%%% FUNCTION DEFINITIONS BELOW THIS LINE

%% This auxilliary function calculates the fit dataset
function [res,Dfit,CConc,Sfit] = FitFunc(p,t,Kmodel,C0,Dexp)
    % Read parameters off the parameter vector (p)
    t0      = p(1);
    FWHM    = p(2);
    
    % Build the K matrix. Note that in the fit we are using k=1/tau
    Kmat    = Kmodel(1./p(3:end));
    
    
    IRFfnc  = @(T,T0,W) (2./W)*sqrt(log(2)/pi)*exp(-4*log(2).*((T-T0)./W).^2)';
    Conc_fn = @(t) kineticsKmat_simu(t,Kmat,C0,0,[]);
    
    % Use the modified IRFconvol routine (IRFconvol_SP), which sets the concentration of the 1st species to 1 at t<t0
    % To use the standard routine, replace by IRFconvol. See manuscript and SI for details.
    CConc   = IRFconvol_SP(IRFfnc,t,t0,FWHM,Conc_fn); 
    
    % Calculate Sfit and Dfit matrices (linear least squares)
    Sfit    = CConc\Dexp;
    Dfit    = CConc*Sfit;
    
    % Calculate sum of squares of residuals (SSR) and non-negative SAS penalty function (NNR)
    % Set the pre-factor of the NNR to 0 to disable it. See manuscript and SI for details.
    NNR     = sum(Sfit(Sfit<0));
    % Tuning parameter for penalised least-squares fit
    GAMMA   = 0; 
    % Penalised norm of residuals (see SI for discussion)
    res     = norm(Dexp-Dfit).^2 + GAMMA*NNR.^2;

end

%% This function generates a dark blue/white/red colour map for visualisation of positive/negative contour plots
function newmap = darkb2r(cmin_input,cmax_input,n_total,n_whites)
    % check the input
    if nargin > 4
       warning('Input two variables, the range of caxis , for example : colormap(darkb2r(-3,3)) or three, darkb2r(-3,3,ncolors)');
    elseif nargin > 2 && nargin < 4
       n_total = 250;
       n_whites = 1;
    end
    
    if cmin_input >= cmax_input
        error('Input error: The color range must be from a smaller one to a larger one');
    end
    
    %% control the figure caxis 
    lims = [0 1]; % Get figure caxis formation
    
    %% color configuration : from blue to to white then to red
    red_top     = [0.5 0 0];
    red_middle  = [1 0 0];
    white_middle= [1 1 1];
    blue_middle = [0 0 1];
    blue_bottom = [0 0 0.5];
    
    if mod(n_whites,2) ~= 0
        n_whites    = n_whites - 1;
    end
    
    n_reds      = round((n_total - n_whites)/2)+1;
    n_blues     = round((n_total - n_whites)/2)+1;
    
    %% Color interpolation 
    m=2/3; % Position of the set colors
    color_input     = [blue_bottom; blue_middle; white_middle];
    oldsteps        = [-1 -m 0];
    newsteps        = linspace(-1, 0, n_blues);  
    bluemap         = zeros(n_blues,3);
    for j=1:3
       bluemap(:,j) = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
    end
    
    color_input     = [white_middle; red_middle; red_top];
    oldsteps        = [0 m 1];
    newsteps        = linspace(0, 1, n_reds);  
    redmap          = zeros(n_reds,3);
    for j=1:3
       redmap(:,j)  = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
    end
    
    bluemap(end,:)  = [];
    redmap(1,:)     = [];
    
    if n_whites ~= 0
        whitemap        = ones(n_whites,3);
        newmap_all      = [bluemap;whitemap;redmap];
    else
        newmap_all      = [bluemap;redmap];
    end
    
    if (cmin_input < 0)  &&  (cmax_input > 0)  
            
        if abs(cmin_input) < cmax_input 
             
            % |--------|---------|--------------------|    
          % -cmax      cmin       0                  cmax         [cmin,cmax]
     
           start_point = max(round((cmin_input+cmax_input)/2/cmax_input*n_total),1);
           newmap = squeeze(newmap_all(start_point:n_total,:));
           
        elseif abs(cmin_input) >= cmax_input
            
             % |------------------|------|--------------|    
           %  cmin                0     cmax          -cmin         [cmin,cmax]   
           
           end_point = max(round((cmax_input-cmin_input)/2/abs(cmin_input)*n_total),1);
           newmap = squeeze(newmap_all(1:end_point,:));
        end
           
    elseif cmin_input >= 0
    
           if lims(1) < 0 
               disp('caution:')
               disp('there are still values smaller than 0, but cmin is larger than 0.')
               disp('some area will be in red color while it should be in blue color')
           end
            % |-----------------|-------|-------------|    
          % -cmax               0      cmin          cmax         [cmin,cmax]
     
           start_point = max(round((cmin_input+cmax_input)/2/cmax_input*n_total),1);
           newmap = squeeze(newmap_all(start_point:n_total,:));
    
    elseif cmax_input <= 0
           if lims(2) > 0 
               disp('caution:')
               disp('there are still values larger than 0, but cmax is smaller than 0.')
               disp('some area will be in blue color while it should be in red color')
           end
           
             % |------------|------|--------------------|    
           %  cmin         cmax    0                  -cmin         [cmin,cmax]      
           end_point = max(round((cmax_input-cmin_input)/2/abs(cmin_input)*n_total),1);
           newmap = squeeze(newmap_all(1:end_point,:));
    end
end