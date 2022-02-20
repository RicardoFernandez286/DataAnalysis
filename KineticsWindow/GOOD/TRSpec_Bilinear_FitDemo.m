% Description:  This script fits time-dependent spectra based on coupled first-order reactions.
%               Useful for testing the robustness of a given model to noise.
%
% Tested and implemented in MATLAB R2021b
% v1.0

% Disable warnings and create wait bar to show progress
warning('off','MATLAB:rankDeficientMatrix');    % Disable warnings about rank-defficient matrices.
warning('off','MATLAB:illConditionedMatrix');   % Disable warnings about near-singular matrices.
wb = waitbar(0,'Fitting simulated data...');

% Define simulation parameters
% (Kinetic parameters are defined inside TRSpec_Bilinear_Simu)
NoisePercent= 2;
Nfits       = 50;

% Simulate data once to get array sizes
[WL,t,Dexp,St] = TRSpec_Bilinear_Simu(0,NoisePercent);

kTrue = [0.5 0.3 0.4 0.1]; % True values of the rate constants

%% Initialise fit parameters and define model
C0 = [1 0 0 0]';
p0 = 1./[0.5 0.3 0.4 0.1];
% Set the upper and lower bounds for the fitted rates
LB = 1e-3*ones(1,length(p0));
UB = 2000*ones(1,length(p0));

%%% Define the model here, by editing the shape of the K matrix
Kmodel   = @(k) [-k(1)  0       0           0 
                k(1)    -k(2)   k(3)        0
                0       k(2)    -k(4)-k(3)  0
                0       0       k(4)        0];  

nC  = size(St,1); % Number of species
nT  = length(t);
nWL = length(WL);

% Preallocate arrays
pFit        = zeros(Nfits,nC);
SASfit      = zeros(nC,nWL,Nfits);
Cfit        = zeros(nT,nC,Nfits);
chi2        = zeros(Nfits,1);

tic;
for i=1:Nfits  
    %% Simulate Data
    [WL,t,Dexp,St] = TRSpec_Bilinear_Simu(0,NoisePercent);
    
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
                        'Display','off');

    [pFit(i,:),resnorm,~,~,output,~,J] = lsqnonlin(@(p) FitFunc(p,t,Kmodel,C0,Dexp),p0,LB,UB,options);
    [~,Dbest,Cfit(:,:,i),SASfit(:,:,i)] = FitFunc(pFit,t,Kmodel,C0,Dexp);
    res     = Dexp-Dbest;
    chi2(i) = resnorm./(numel(res)-length(pFit)-1);
    waitbar(i./Nfits,wb,['Fitting dataset ' num2str(i) ' of ' num2str(Nfits)])
end
delete(wb);
toc
%%
% Plot results
fh=figure(1);
clf(fh);
fh.Color = 'w';
ax = axes('parent',fh);

KtruePlot = kTrue.*ones(size(t));
cmap = turbo(nC+1);

hold(ax,'on')
for j=1:nC
    plot(ax,1:Nfits,1./pFit(:,j),'LineWidth',2,'Color',cmap(j,:))
    plot(ax,KtruePlot(:,j),'--','LineWidth',1.5,'Color',cmap(j,:),'HandleVisibility','off')
end
hold(ax,'off')
xlabel(ax,'Fit Number','FontWeight','bold');
ylabel(ax,'Rate constant (s^{-1})','FontWeight','bold');
legend(ax,{'k_{1}','k_{2f}','k_{2b}','k_{3}'},'Location','best','orientation','horizontal')
ax.FontSize =16;
xlim(ax,[1,Nfits]);
box(ax,'on');

%%%%%% END OF FILE

%%%%%% FUNCTION DEFINITIONS BELOW THIS LINE
%% This auxilliary function calculates the fit dataset
function [res,Dfit,CConc,Sfit] = FitFunc(p,t,Kmodel,C0,Dexp)
% Build the K matrix. Note that in the fit we are using k=1/tau
Kmat    = Kmodel(1./p);

Conc_fn = @(t) kineticsKmat_simu(t,Kmat,C0,0,[]);
CConc   = Conc_fn(t); 

% Calculate Sfit and Dfit matrices (linear least squares)
Sfit    = CConc\Dexp;
Dfit    = CConc*Sfit;

% Calculate sum of squares of residuals (SSR) and non-negative SAS penalty function (NNR)
% Set the pre-factor of the NNR to 0 to disable it. See manuscript and SI for details.
NNR     = sum(Sfit(Sfit<0));
res     = norm(Dexp-Dfit).^2 + 0*NNR.^2;

end