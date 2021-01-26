% scriptdir   = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\20) M-H solvation - na\Data\CLS';
% subdir      = 'IrHCOP3';

scriptdir = 'E:\Masterarbeit\Data\Lab 4\IrHP3\CLS';
subdir = [];

filelist    = dir([scriptdir filesep subdir]);

doFit = 1;

%% Parse the names into concentrations and make a list
names       = {filelist.name}';
names       = flipud(names(contains(names,'.csv')));
% names       = names(1:end-1);
Nsolvents   = length(names);
solventPos  = 3:5;

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
    theta(i,1)      = solvents.KT(svt_ID,3);
    solventNames(i) = solvents.Names(svt_ID);
end

[theta,idx]     = sort(theta,'ascend');
names           = names(idx);
solventNames    = solventNames(idx);

%% Read the data
for i=1:Nsolvents
    data        = readmatrix([scriptdir filesep subdir filesep names{i}]);

    t2delays{i} = data(:,1);
    Ncols       = size(data,2)-1;

    SpecDiff{i} = data(:,2:end);
    
%     t2plot = 10;
%     id = findClosestId2Val(t2delays{i},t2plot);
%     pepita(i,:) = SpecDiff{i}(id,:);
%     pepito(i,:) = SpecDiff{i}(1,:);
end

%% Fit options
fitfun = @(P,t)P(1).*exp(-t./P(2)) + P(3);

start_param = [0.8 3 0.1; 0.6 3 0.2];
UB_all      = [1   10 1;  1   10  1];
LB_all      = [0  0.25 0; 0  0.25 0];  

t2start_fit = [0.25 0.25]; 
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
            
            [fitted_param(i,:,q),resnorm,residuals,exitflag,output,lambda,jacobian_fit] = lsqcurvefit(fitfun,fitPar0,xdata,ydata,LB,UB,options);
            ci = nlparci(fitted_param(i,:,q),residuals,'jacobian',jacobian_fit);
            fitted_err(i,:,q) = ci(:,2) - ci(:,1);
            
            yFit = fitfun(fitted_param(i,:,q),xdata);
            plot(ax,xdata,yFit,'-','Color',cmap(i,:),'DisplayName',solventNames(i),'LineWidth',2);
        else
            plot(ax,t2delays{i},SpecDiff{i}(:,q),'o-','Color',cmap(i,:),'DisplayName',solventNames(i),'LineWidth',2);
        end
    end
    hold(ax,'off');
    
    ylim(ax,[0 1.1]);
    yline(ax,1,'handlevisibility','off');
    box(ax,'on');
    legend(ax,'show');
end
