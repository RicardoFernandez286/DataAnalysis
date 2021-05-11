clear all;
% scriptdir   = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\20) M-H solvation - na\Data\CLS';
% subdir      = 'IrHCOP3';

% scriptdir = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\20) M-H solvation - na\Data\CLS';
% subdir = 'IrHCOP3_old';

rootdir = 'E:\Masterarbeit\Data\Lab 4\VC-H2 high frequency';
% rootdir = 'E:\Masterarbeit\Data\Lab 4\VC-H2';
% rootdir = 'E:\Masterarbeit\Data\Lab 4\IrHP3';

figdir  = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\20) M-H solvation - na\Manuscript\Figures\MATLAB';
complexname = 'VCH2';

subdir  = 'CLS';

filelist = dir([rootdir filesep subdir]);

doFit = 1;
errbar = 0;

PlotScale=3; % 3=theta, 4=ET30, 5=Viscosity
% parplot = 3;  % 1 = A1, 2 = tau, 3 = A0


%% Parse the names into concentrations and make a list
names       = {filelist.name}';
names       = flipud(names(contains(names,'.csv')));
% names       = names(1:end-1);
Nsolvents   = length(names);
solventPos  = 2:5;

% Initialise variables
Solvent         = strings(Nsolvents,1);
SpecDiff        = cell(Nsolvents,1);
t2delays        = cell(Nsolvents,1);

% Read the solvents and sort them by polarity (Theta scale)
load('Solvents')

fn_parts = cell(Nsolvents,1);
theta    = zeros(Nsolvents,1);
solventNames = string;

for i=1:Nsolvents
    fn_parts{i}     = strsplit(names{i},'_');
    svt_ID          = matches(solvents.Names(:,2),fn_parts{i}(solventPos),'IgnoreCase',true);
    theta(i,1)      = solvents.KT(svt_ID,PlotScale);
    solventNames(i) = solvents.Names(svt_ID);
end

[theta,idx]     = sort(theta,'ascend');
names           = names(idx);
solventNames    = solventNames(idx);

%% Read the data
for i=1:Nsolvents
    data        = readmatrix([rootdir filesep subdir filesep names{i}]);

    t2delays{i} = data(:,1);
    Ncols       = size(data,2)-1;

%     Amix        = [0.81 0.19; 0.3 0.7];
%     SpecDiff{i} = (Amix\data(:,2:end)')';
   
    SpecDiff{i} = data(:,2:end);
    
%     t2plot = 10;
%     id = findClosestId2Val(t2delays{i},t2plot);
%     pepita(i,:) = SpecDiff{i}(id,:);
%     pepito(i,:) = SpecDiff{i}(1,:);
end

%% Fit options
fitfun = @(P,t)P(1).*exp(-t./P(2)) + P(3);

start_param = [0.8 3 0.1; 0.6 3 0.2];
UB_all      = [1   15 1;  1   15  1];
LB_all      = [0  0.25 0; 0  0.25 0];  

t2start_fit = [0. 0.]; 
t2end_fit   = [25 15];

%% Plot the data
cmap = othercolor('Mrainbow',Nsolvents);

for q=1:Ncols
    % Create figure
    fh          = figure(q); clf(fh);
    fh.Units    = 'normalized';
    ax          = axes('parent',fh);

    hold(ax,'on');
    for i=1:Nsolvents    
        if doFit == 1
            id_s  = findClosestId2Val(t2delays{i},t2start_fit(q));
            id_e  = findClosestId2Val(t2delays{i},t2end_fit(q));
            xdata = t2delays{i}(id_s:id_e);
            ydata = SpecDiff{i}(id_s:id_e,q);
            
            fitPar0 = start_param(q,:);
            LB      = LB_all(q,:);
            UB      = UB_all(q,:);
            
            options = optimoptions('lsqcurvefit',...
                'MaxFunctionEvaluations',1000,...
                'MaxIterations',1000,...
                'Algorithm','trust-region-reflective',...  %'levenberg-marquardt' 'trust-region-reflective'
                'OptimalityTolerance',1e-15,...
                'FunctionTolerance',1e-15,...
                'TypicalX',fitPar0,...
                'SubproblemAlgorithm','factorization',...
                'StepTolerance',1e-15);%,...
    
            plot(ax,xdata,ydata,'o','Color',cmap(i,:),'DisplayName',solventNames(i),'LineWidth',2,'handlevisibility','off');
            
            [Fit_Par(i,:,q),resnorm,residuals,exitflag,output,lambda,jacobian_fit] = lsqcurvefit(fitfun,fitPar0,xdata,ydata,LB,UB,options);
            ci = nlparci(Fit_Par(i,:,q),residuals,'jacobian',jacobian_fit);
            Fit_Err(i,:,q) = ci(:,2) - ci(:,1);
            
            xinterp = linspace(min(xdata),max(xdata),1000);
            yFit = fitfun(Fit_Par(i,:,q),xinterp);
            plot(ax,xinterp,yFit,'-','Color',cmap(i,:),'DisplayName',solventNames(i),'LineWidth',2);
            
            title(ax,['\nu_{' num2str(q+1) '}']);
        else
            plot(ax,t2delays{i},SpecDiff{i}(:,q),'o-','Color',cmap(i,:),'DisplayName',solventNames(i),'LineWidth',2);
        end
    end
    hold(ax,'off');
    
    ylim(ax,[0 1.1]);
    yline(ax,1,'handlevisibility','off');
    box(ax,'on');
    ax.FontSize=18;
    legend(ax,'show');
%     legend(ax,'location','best');
    legend(ax,'NumColumns',3);
    legend(ax,'fontsize',12);
    ylabel(ax,'Centerline Slope (CLS)','FontWeight','bold');
    xlabel(ax,'Population Delay (t_{2}, ps)','FontWeight','bold');
    xlim(ax,[0,25]);
end
%% Resize figures
fh1 = figure(1);
fh2 = figure(2);
for fig=[fh1,fh2]
    fig.Units = 'pixels';
    fig.Color = 'w';
end
fh1.Position(3:4) = [640 400];

for fig=[fh1,fh2]
    fig.Position(3:4) = fh1.Position(3:4);
end

%% Plot the fit parameters across solvents
for p=1:3
    parplot=p;
    % Parse axes labels
    switch parplot
        case 1
            y_string = 'A_{1} / CLS Amplitude';
            fig_suf = 'A1';
        case 2
            y_string = '\tau_{C} (ps)';
            fig_suf = 'tauC';
        case 3
            y_string = 'A_{0} / CLS Offset';
            fig_suf = 'A0';
    end

    switch PlotScale
        case 3
            x_string = 'Solvent Acidity/Basicity, \theta = \alpha-\beta';
            fig_pref = 'theta';
        case 4
            x_string = 'E_{T}(30) Polarity';
            fig_pref = 'ET30';
        case 5
            x_string = '\eta (mPa s)';
            fig_pref = 'visc';
    end

    fh3(p) = figure(2+p);
    fh3(p).Color = 'w';
    clf(fh3(p));
    ax2 = axes('parent',fh3(p));

    hold(ax2,'on');
    if errbar == 1
        errorbar(ax2,theta,squeeze(Fit_Par(:,parplot,1))',squeeze(Fit_Err(:,parplot,1))','o-','Linewidth',1.5,'Color','b','DisplayName','\nu_{2}')
        errorbar(ax2,theta,squeeze(Fit_Par(:,parplot,2))',squeeze(Fit_Err(:,parplot,2))','o-','Linewidth',1.5,'Color','r','DisplayName','\nu_{3}')
    else
        plot(ax2,theta,squeeze(Fit_Par(:,parplot,1))','sq-','Linewidth',1.5,'Color','b','MarkerSize',8,'DisplayName','\nu_{2}')
        plot(ax2,theta,squeeze(Fit_Par(:,parplot,2))','o -','Linewidth',1.5,'Color','r','MarkerSize',8,'DisplayName','\nu_{3}')
    end

    hold(ax2,'off');
    legend(ax2,'show');
    legend(ax2,'location','best');
    % legend(ax2,'orientation','horizontal');
    box(ax2,'on');
    xlabel(ax2,x_string,'FontWeight','bold');
    ylabel(ax2,y_string,'FontWeight','bold');
    ax2.FontSize = 18;
    axis(ax2,'tight');

    if PlotScale==3
        xline(ax2,0,'Handlevisibility','off');
    end
    title(ax2,fig_suf);
    
    fh3(p).Position(3:4) = fh1.Position(3:4);
    fh3(p).Units = 'pixels';
    fh3(p).Color = 'w';
    fh3(p).Position(3:4) = fh1.Position(3:4);
    fh3(p).Units = 'normalized';
    
    exportgraphics(fh3(p),[figdir filesep complexname '_CLS_' fig_pref '_' fig_suf '.eps'])
end

%% Export
[theta,squeeze(Fit_Par(:,parplot,1)),squeeze(Fit_Par(:,parplot,2))]