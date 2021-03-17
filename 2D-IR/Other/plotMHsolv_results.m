clear all;

DoFit   = 1;
DoSave  = 0;

NormType='Sqrt'; % 'Old' 'Sqrt' 'Prod' 'Squared'

WhatPlot='IrHCOP3';

%% Get a list of MAT files
switch WhatPlot
    case 'IrHCOP3'
        scriptdir   = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\20) M-H solvation - na\Data\FitResults';
%         scriptdir   = 'E:\Masterarbeit\Data\Lab 4\IrHP3\CLS';
        subdir      = 'IrHCOP3';
end

% plotWhat    = 'Xpeak ESA'; % Xpeak or Diagonal + GSB/ESA
plotWhat    = 'Diagonal ESA'; % Xpeak or Diagonal + GSB/ESA
% plotWhat    = 'SpecDiff';
plotFormat  = 'Vertical'; % 'Horizontal' or 'Vertical'

xpos        = 0.83; % 0.8 for Re18; 0.83 for CNBz

% xpos        = 1.04; % 1.04 for Re18; 1.06 for CNBz; 0.1 horizontal
% xpos        = -0.175;
ypos        = 0.35;

FontSize    = 14;
LineWidth   = 1;
MarkerSize  = 5;
line_up     = 'o';
line_dw     = 'o';
    
XScale      = 'lin';
YmaxPlot    = 'p';

filelist    = dir([scriptdir filesep subdir]);

%% Parse the names into concentrations and make a list
names       = {filelist.name}';
names       = flipud(names(contains(names,'.mat')));
% names       = names(1:end-1);
Nconc       = length(names);
solventPos  = 3;

% Initialise variables
Solvent         = strings(Nconc,1);
VolumeData_GSB  = cell(Nconc,1);
VolumeData_ESA  = cell(Nconc,1);
SpecDiff        = cell(Nconc,1);
SSR_fit         = zeros(Nconc,1);
StepSize        = zeros(Nconc,1);
delays          = cell(Nconc,1);
xpeak_UP        = cell(Nconc,1);
xpeak_DW        = cell(Nconc,1);

% Create figure
fh              = figure(1);
clf(fh);
fh.Units        = 'normalized';
% fh.Position(2)  = 0.1;
% fh.Position(4)  = fh.Position(4)*2;
switch plotFormat
    case 'Vertical'
%         fh.Position     = [0.3536    0.1000    0.275    0.60];
        ax_UP           = subplot(2,1,1);
        ax_DW           = subplot(2,1,2);
    case 'Horizontal'
%         fh.Position     = [0.15    0.15    0.6    0.45];
        ax_UP           = subplot(1,2,1);
        ax_DW           = subplot(1,2,2);
end

fh.Color        = [1 1 1];
fh.Units        = 'pixels';
ax_UP.FontSize  = FontSize;
ax_DW.FontSize  = FontSize;

box(ax_UP,'on');
box(ax_DW,'on');

hold(ax_UP,'on');
hold(ax_DW,'on');

cmFW = colormap(ax_UP,(othercolor('Mrainbow',Nconc)));
cmBW = colormap(ax_DW,(othercolor('Mrainbow',Nconc)));

% cmFW = colormap(ax_UP,flip(othercolor('Blues3',Nconc)));
% cmBW = colormap(ax_DW,flip(othercolor('Reds3',Nconc)));

%% Read the MAT files one by one and store the results in a cell array
% Read the solvents and sort them by polarity (Theta scale)
load('Solvents')

fn_parts = cell(Nconc,1);
theta    = zeros(Nconc,1);
solventNames = string;

for i=1:Nconc
    fn_parts{i} = strsplit(names{i},'_');
    svt_ID = matches(solvents.Names(:,2),fn_parts{i}{solventPos},'IgnoreCase',true);
    theta(i,1) = solvents.KT(svt_ID,3);
    solventNames(i) = solvents.Names(svt_ID);
end
% svt = 

[theta,idx]     = sort(theta,'ascend');
names           = names(idx);

% Now plot the stuff
for i=1:Nconc
    load([scriptdir filesep subdir filesep names{i}])

    %%% REDO THE NORMALISATION.
        % The volume structure dimensions are [time x peak # x (GSB/ESA)]
    switch NormType
        case 'Sqrt'
            % Normalise by dividing by sqrt(pump*probe)
            for p=1:2 % 1=GSB, 2=ESA ---- VOLUMES
                NormVols(:,1,p) = ((-1)^p)*fitPar.Vols(:,1,p)./max(abs(fitPar.Vols(:,1,p))); % 13CO diag
                NormVols(:,2,p) = ((-1)^p)*fitPar.Vols(:,2,p)./max(abs(fitPar.Vols(:,2,p))); % 12CO diag
                NormVols(:,3,p) = ((-1)^p)*fitPar.Vols(:,3,p)./sqrt(max(abs(fitPar.Vols(:,1,p))).*max(abs(fitPar.Vols(:,2,p)))); % 13->12 uphill Xpeak
                NormVols(:,4,p) = ((-1)^p)*fitPar.Vols(:,4,p)./sqrt(max(abs(fitPar.Vols(:,1,p))).*max(abs(fitPar.Vols(:,2,p)))); % 12->13 downhill Xpeak
            end
           
            for p=1:2 % 1=GSB, 2=ESA  ---- ERRORS
                NormErr(:,1,p)  = ((-1)^p)*fitErr.Vols(:,1,p)./max(abs(fitErr.Vols(:,1,p))); % 13CO diag
                NormErr(:,2,p)  = ((-1)^p)*fitErr.Vols(:,2,p)./max(abs(fitErr.Vols(:,2,p))); % 12CO diag
                NormErr(:,3,p)  = ((-1)^p)*fitErr.Vols(:,3,p)./sqrt(max(abs(fitErr.Vols(:,1,p))).*max(abs(fitErr.Vols(:,2,p)))); % 13->12 uphill Xpeak
                NormErr(:,4,p)  = ((-1)^p)*fitErr.Vols(:,4,p)./sqrt(max(abs(fitErr.Vols(:,1,p))).*max(abs(fitErr.Vols(:,2,p)))); % 12->13 downhill Xpeak
            end
    end   
    VolumeData_GSB{i}   = NormVols(:,:,1);
    VolumeData_ESA{i}   = NormVols(:,:,2);
    SpecDiff{i}         = fitPar.C;
    
    X0Pos(:,i)      = fitPar.X0;
    Anharm(:,i)     = fitPar.Anharm;
    AMPRATIO(i,:)   = squeeze(fitPar.Amps(1,1:2,2).*sqrt(fitPar.Sx(1:2))'.*sqrt(fitPar.Sy(1:2))');
    AmpQuot = AMPRATIO(:,2)./AMPRATIO(:,1);
    
 %%   
    delays{i}     = t2delays;
    
    decay = ones(size(t2delays));
    
    xpeak_UP{i}   = NormVols(:,3,2).*decay;
    xpeak_DW{i}   = NormVols(:,4,2).*decay;
    
    SSR_fit(i)       = SSR;
    StepSize(i)      = output_st.stepsize;
    switch plotWhat
        case 'Diagonal GSB'
            plotData_up = NormVols(:,1,1).*decay;
            plotData_dw = NormVols(:,2,1).*decay;
        case 'Xpeak GSB'
            plotData_up = NormVols(:,3,1).*decay;
            plotData_dw = NormVols(:,4,1).*decay;
        case 'Diagonal ESA'
            plotData_up = NormVols(:,1,2).*decay;
            plotData_dw = NormVols(:,2,2).*decay;
        case 'Xpeak ESA'
            plotData_up = NormVols(:,3,2).*decay;
            plotData_dw = NormVols(:,4,2).*decay;
        case 'SpecDiff'
            plotData_up = fitPar.C(:,1);
            plotData_dw = fitPar.C(:,2);
        case 'Xpeak Diff'
            plotData_up = (NormVols(:,3,1)+NormVols(:,3,2)).*decay;
            plotData_dw = (NormVols(:,4,1)+NormVols(:,4,2)).*decay;
    end
    
    plot(ax_UP,t2delays,plotData_up,line_up,'MarkerSize',MarkerSize,'LineWidth',LineWidth,'Color',cmFW(i,:),'HandleVisibility','off');
    plot(ax_DW,t2delays,plotData_dw,line_dw,'MarkerSize',MarkerSize,'LineWidth',LineWidth,'Color',cmBW(i,:),'HandleVisibility','off');
    
    plot(ax_UP,NaN,NaN,'-','LineWidth',4,'Color',cmFW(i,:),'DisplayName',solventNames(i));
    switch plotFormat
    case 'Vertical'
        plot(ax_DW,NaN,NaN,'-','LineWidth',4,'Color',cmBW(i,:),'DisplayName',solventNames(i))
    end
end
% Customize axes
xpeakUp     = 'Re(^{13}CO) \rightarrow Re(^{12}CO) \rm{\times10}';
xpeakDown   = 'Re(^{13}CO) \leftarrow Re(^{12}CO) \rm{\times10}';
titUP       = ['Uphill energy transfer' ', ' xpeakDown];
titDW       = ['Downhill energy transfer' ', ' xpeakUp];

% titUP       = 'Uphill energy transfer';
% titDW       = 'Downhill energy transfer';

% Set axes limits
axis(ax_UP,'tight');
axis(ax_DW,'tight');

%%% Nice formatting
% Titles
title(ax_UP,titUP,'FontSize',FontSize);
title(ax_DW,titDW,'FontSize',FontSize);

% Axis labels
xlabel(ax_UP,'t_2 delay (ps)','FontWeight','bold','FontSize',FontSize);
xlabel(ax_DW,'t_2 delay (ps)','FontWeight','bold','FontSize',FontSize);

if contains(plotWhat,'Xpeak') % Xpeaks
    ylim(ax_UP,[-0.05 0.8]);
    ylim(ax_DW,[-0.05 0.4]);
    ylabel(ax_UP,'Normalised Volume','FontWeight','bold','FontSize',FontSize);
%     ylabel(ax_UP,'Population (%)','FontWeight','bold','FontSize',FontSize);
    switch plotFormat
        case 'Vertical'
            ylabel(ax_DW,'Normalised Volume','FontWeight','bold','FontSize',FontSize);
%             ylabel(ax_DW,'Population (%)','FontWeight','bold','FontSize',FontSize);
    end
elseif contains(plotWhat,'SpecDiff') % Spectral diffusion
    ylim(ax_UP,[-0.05 1]);
    ylim(ax_DW,[-0.05 1]);
    ylabel(ax_UP,'C (t_{2})','FontWeight','bold','FontSize',FontSize);
    switch plotFormat
        case 'Vertical'
            ylabel(ax_DW,'C (t_{2})','FontWeight','bold','FontSize',FontSize);
    end
else % Diagonal
    ylim(ax_UP,[-0.05 1]);
    ylim(ax_DW,[-0.05 1]);
    ylabel(ax_UP,'Normalised Volume','FontWeight','bold','FontSize',FontSize);
%     ylabel(ax_UP,'Norm. Vol.','FontWeight','bold','FontSize',FontSize);
    switch plotFormat
        case 'Vertical'
            ylabel(ax_DW,'Normalised Volume','FontWeight','bold','FontSize',FontSize);
%             ylabel(ax_DW,'Population (%)','FontWeight','bold','FontSize',FontSize);
    end
end


% Axis limits
if isa(YmaxPlot,'double')
    ylim(ax_UP,[0 YmaxPlot]);
    ylim(ax_DW,[0 YmaxPlot]);
else
    axis(ax_UP,'tight');
    axis(ax_DW,'tight');
end

xlim(ax_UP,[0 40]);
xlim(ax_DW,[0 40]);


% % Add zero line
% hline_FW = yline(ax_UP,0,'HandleVisibility','off'); hline_FW.Color = [0.5 0.5 0.5];
% hline_BW = yline(ax_DW,0,'HandleVisibility','off'); hline_BW.Color = [0.5 0.5 0.5];

% Legends
switch plotFormat
    case 'Vertical'
        plot(ax_UP,NaN,NaN,'.','Color','w','DisplayName','')
%         text(ax_DW,xpos,ypos,['\bf{% ' diluent '}'],'FontSize',12,'FontWeight','bold','Units','normalized')
        plot(ax_DW,NaN,NaN,'.','Color','w','DisplayName','')
        leg_DW = legend(ax_DW,'show');
        legend(ax_DW,'boxoff','FontWeight','bold')
        legend(ax_DW,'location','eastoutside')
        leg_DW.ItemTokenSize = [20,10];
        leg_UP = legend(ax_UP,'show');
        legend(ax_UP,'boxoff','FontWeight','bold')
        legend(ax_UP,'location','bestoutside')
        leg_UP.ItemTokenSize = [20,10];
        delete(leg_DW);
    case 'Horizontal'
        text(ax_UP,xpos,1.2,['\bf{% ' diluent '}'],'FontSize',12,'FontWeight','bold','Units','normalized')
        leg_UP = legend(ax_UP,'show');
        legend(ax_UP,'boxoff')
        legend(ax_UP,'FontWeight','bold')
        legend(ax_UP,'orientation','horizontal')
        leg_UP.Position = [0.3 0.9 0.5 0.06];
        leg_UP.ItemTokenSize = [20,10];
end

% Make figure resizable
fh.Units        = 'normalized';

switch plotFormat
    case 'Vertical'
        ax_UP.Position = [0.10 0.58 0.70  0.35];
        ax_DW.Position = [0.10 0.10 0.70  0.35];
        fh.Position(3:4)    = [0.240  0.51];
        xlabel(ax_UP,[]);
        title(ax_UP,'(A) Upper Cross Peak','Units','normalized','position',[0.3 1.05]);
        title(ax_DW,'(B) Lower Cross Peak','Units','normalized','position',[0.33 1.05]);
%         title(ax_UP,'(A) ^{13}CO diagonal peak','Units','normalized','position',[0.25 1.05]);
%         title(ax_DW,'(B) ^{12}CO diagonal peak','Units','normalized','position',[0.25 1.05]);
        leg_UP.Position(1:2) = [0.85  (0.5-leg_UP.Position(4)/2)];
        annotation(fh,'textbox',[xpos -0.05+leg_UP.Position(4)+leg_UP.Position(2) 0.5 0.05],'String',['\bf{Solvent:}'],'FontSize',12,'FontWeight','bold','Units','normalized','EdgeColor','none');
    case 'Horizontal'
        ax_UP.Position = [0.10 0.15 0.385 0.65];
        ax_DW.Position = [0.55 0.15 0.385 0.65];
end
hold(ax_UP,'off');
hold(ax_DW,'off');

ax_UP.XScale = XScale;
ax_DW.XScale = XScale;

% ax_UP.XTick = [0.1 1 10 100 1000 10000];
% ax_DW.XTick = [0.1 1 10 100 1000 10000];

% ax_UP.XTick = [0.1 1 10 100 1000 10000];
% ax_DW.XTick = [0.1 1 10 100 1000 10000];


ax_UP.TickLength = 2*ax_UP.TickLength;
ax_DW.TickLength = 2*ax_DW.TickLength;
leg_UP.Position(2) = leg_UP.Position(2)-0.04;

% FIGURE SIZE
fh.Units = 'pixels';
% fh.Position = [275 718 1285 620];
fh.Position(3) = 614.4*3/4;

%%%
% DoFit = 0;
if DoFit == 0
 return
end

%% Do kinetic model fit
Xpeaks  = contains(plotWhat,'Xpeak','ignorecase',true);
C_t     = contains(plotWhat,'SpecDiff','ignorecase',true);

if Xpeaks == 1
    fitPar0 = [-1  0.25     20      0.01];
    UB      = [1   15       100     1];
    LB      = [-1  0.25     0.25   	-1];
    Nexp    = 2.1;
elseif C_t == 1
    fitPar0 = [1   0.25     0.01];
    UB      = [1   15       1];
    LB      = [0   0.25     0];
    Nexp    = 1.1;
    
    CLS_lastDelay   = 20;
    lastDelay_ID    = findClosestId2Val(t2delays,CLS_lastDelay);
    t2delays        = t2delays(1:lastDelay_ID);
else
    fitPar0 = [1   0.25];
    UB      = [1   15  ];
    LB      = [0   0.25];
    Nexp    = 1;
end

lastDelay_ID = length(t2delays);


for i=1:Nconc
    if contains(plotWhat,'SpecDiff','ignorecase',true)
        [Fit_Par(:,:,i),Fit_Curves(:,:,i),Fit_Err(:,:,i)] = FitDoubleExp(t2delays,SpecDiff{i}(1:lastDelay_ID,:),fitPar0,LB,UB,Nexp);
    else
        [Fit_Par(:,:,i),Fit_Curves(:,:,i),Fit_Err(:,:,i)] = FitDoubleExp(t2delays,VolumeData_ESA{i}(1:lastDelay_ID,:),fitPar0,LB,UB,Nexp);
    end
end


hold(ax_UP,'on');
hold(ax_DW,'on');

for i=1:Nconc
    plot(ax_UP,t2delays,Fit_Curves(:,1,i),'-','Color',cmFW(i,:),'LineWidth',LineWidth,'HandleVisibility','off');
    plot(ax_DW,t2delays,Fit_Curves(:,2,i),'-','Color',cmBW(i,:),'LineWidth',LineWidth,'HandleVisibility','off');
end

hold(ax_UP,'off');
hold(ax_DW,'off');

fh2 = figure(2);
clf(fh2);
ax2 = axes('parent',fh2);

if C_t == 0
    hold(ax2,'on');
        errorbar(ax2,theta,squeeze(Fit_Par(1,2,:))',squeeze(Fit_Err(1,2,:))','o','Color','b')
        errorbar(ax2,theta,squeeze(Fit_Par(2,2,:))',squeeze(Fit_Err(2,2,:))','o','Color','r')
    hold(ax2,'off');

    m = mean(squeeze(Fit_Par(1,2,:)));
    s = std(squeeze(Fit_Par(1,2,:)));
    disp(['Tau1 =' num2str(m) ' ± ' num2str(s) ' ps']);
    yline(ax2,m,'Color','b','LineWidth',2);

    m = mean(squeeze(Fit_Par(2,2,:)));
    s = std(squeeze(Fit_Par(2,2,:)));
    disp(['Tau2 =' num2str(m) ' ± ' num2str(s) ' ps']);
    yline(ax2,m,'Color','r','LineWidth',2);

    if Xpeaks == 1
        m = mean(squeeze(Fit_Par(1,3,:)));
        s = std(squeeze(Fit_Par(1,3,:)));
        disp(['Tau1dec =' num2str(m) ' ± ' num2str(s) ' ps']);

        m = mean(squeeze(Fit_Par(2,3,:)));
        s = std(squeeze(Fit_Par(2,3,:)));
        disp(['Tau2dec =' num2str(m) ' ± ' num2str(s) ' ps']);
    end
else
    m = mean(squeeze(Fit_Par(1,3,:)));
    s = std(squeeze(Fit_Par(1,3,:)));
    disp(['A0_CO =' num2str(m) ' ± ' num2str(s) ' ps']);

    m = mean(squeeze(Fit_Par(2,3,:)));
    s = std(squeeze(Fit_Par(2,3,:)));
    disp(['A0_H =' num2str(m) ' ± ' num2str(s) ' ps']);
    
    hold(ax2,'on');
%         errorbar(ax2,theta,squeeze(Fit_Par(1,3,:))',squeeze(Fit_Err(1,3,:))','o','Color','b')
%         errorbar(ax2,theta,squeeze(Fit_Par(2,3,:))',squeeze(Fit_Err(2,3,:))','o','Color','r')
        plot(ax2,theta,squeeze(Fit_Par(1,3,:))','o-','Color','b','DisplayName','\nu_{1}')
        plot(ax2,theta,squeeze(Fit_Par(2,3,:))','o-','Color','r','DisplayName','\nu_{2}')
    hold(ax2,'off');
    legend(ax2,'show');
    box(ax2,'on');
    xlabel(ax2,'Solvent Acidity/Basicity, \theta = \alpha-\beta','FontWeight','bold');
    ylabel(ax2,'A_{0} / CLS Offset','FontWeight','bold');
    ax2.FontSize = 18;
end
xline(ax2,0,'HandleVisibility','off');
%%%%%%% EOF

function [fitted_param,fitted_curves,fitted_err] = FitDoubleExp(time,VolumeData,fitPar0,LB,UB,Nexp)   
    switch Nexp
        case 1
            fitfun = @(P,t)P(1).*exp(-t./P(2));
            shift  = 0;
        case 1.1
            fitfun = @(P,t)P(1).*exp(-t./P(2)) + P(3);
            shift  = 0;
        case 2
            fitfun = @(P,t)P(1).*exp(-t./P(2)) + (1-P(1)).*exp(-t./P(3)) + P(4);
            shift  = 0;
        case 2.1
            fitfun = @(P,t)P(1).*exp(-t./P(2)) + P(4).*exp(-t./P(3));
            shift  = 2;
    end
    
    options = optimoptions('lsqcurvefit',...
            'MaxFunctionEvaluations',1000,...
            'MaxIterations',1000,...
            'Algorithm','trust-region-reflective',...  %'levenberg-marquardt' 'trust-region-reflective'
            'OptimalityTolerance',1e-15,...
            'FunctionTolerance',1e-15,...
            'TypicalX',fitPar0,...
            'SubproblemAlgorithm','factorization',...
            'StepTolerance',1e-15);%,...
    
    Ncols = 2;
%     Ncols = size(VolumeData,2);
    fitted_param = zeros(Ncols,size(fitPar0,2));
    fitted_err   = zeros(Ncols,size(fitPar0,2));
    
    for j=1:Ncols
        Ydata = VolumeData(:,j+shift);
        [fitted_param(j,:),resnorm,residuals,exitflag,output,lambda,jacobian_fit] = lsqcurvefit(fitfun,fitPar0,time,Ydata,LB,UB,options);
        ci = nlparci(fitted_param(j,:),residuals,'jacobian',jacobian_fit);
        fitted_err(j,:) = ci(:,2) - ci(:,1);
        fitted_curves(:,j) = fitfun(fitted_param(j,:),time);
    end 
end


