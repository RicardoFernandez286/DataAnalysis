%% Get a list of MAT files
scriptdir   = '\\idnetapp-chem.uzh.ch\g_chem_hamm$\Group\Andrea\october_2019';
subdir      = 'fits';

plotWhat    = 'Xpeak GSB + C(t)'; 
% plotWhat    = 'Xpeak ESA'; % Xpeak or Diagonal + GSB/ESA
plotFormat  = 'Horizontal'; % 'Horizontal' or 'Vertical'
concType    = '%'; % '100-%' or '%'
diluent     = 'D_{2}O'; % 'Re(^{13}C^{18}O)'

% diluent     = 'Re(^{13}C^{18}O)';
xpos        = 0; % 1.02 for CNBz, 0.95 for Re18, 0.85 for D2O
                % 0.25 for horizontal

NormType    = 'Sqrt';
DoSave      = 1;

filelist    = dir([scriptdir filesep subdir]);
filelist    = filelist(3:end);

%% Parse the names into concentrations and make a list
names       = {filelist.name}';
Nconc       = length(names);

% Initialise variables
ConcPercent     = zeros(Nconc,1);
PrismID         = strings(Nconc,1);
SolutionID      = strings(Nconc,1);
VoumeData_GSB   = cell(Nconc,1);
VoumeData_ESA   = cell(Nconc,1);
SSR_fit         = zeros(Nconc,1);
StepSize_fit    = zeros(Nconc,1);

% Create figure
fh              = figure(1);
clf(fh);
fh.Units        = 'normalized';
% fh.Position(2)  = 0.1;
% fh.Position(4)  = fh.Position(4)*2;
switch plotFormat
    case 'Vertical'
        fh.Position     = [0.3536    0.1000    0.35    0.75];
        ax_UP           = subplot(2,1,1);
        ax_DW           = subplot(2,1,2);
    case 'Horizontal'
        fh.Position     = [0.15    0.15    0.6    0.45];
        ax_UP           = subplot(1,2,1);
        ax_DW           = subplot(1,2,2);
end

fh.Color        = [1 1 1];
fh.Units        = 'pixels';
ax_UP.FontSize  = 16;
ax_DW.FontSize  = 16;

box(ax_UP,'on');
box(ax_DW,'on');

hold(ax_UP,'on');
hold(ax_DW,'on');

cmUP = colormap(ax_UP,othercolor('Mrainbow',Nconc));
cmDW = colormap(ax_DW,othercolor('Mrainbow',Nconc));

% cmUP = colormap(ax_UP,flip(othercolor('Blues3',Nconc)));
% cmDW = colormap(ax_DW,flip(othercolor('Reds3',Nconc)));

switch plotFormat
    case 'Vertical'
        text(ax_DW,xpos,0.94,['\bf{% ' diluent '}'],'FontSize',12,'FontWeight','bold','Units','normalized')
        plot(ax_DW,NaN,NaN,'.','Color','w','DisplayName','')
        text(ax_UP,xpos,0.94,['\bf{% ' diluent '}'],'FontSize',12,'FontWeight','bold','Units','normalized')
        plot(ax_UP,NaN,NaN,'.','Color','w','DisplayName','')
    case 'Horizontal'
        text(ax_UP,xpos,1.19,['\bf{% ' diluent '}'],'FontSize',12,'FontWeight','bold','Units','normalized')
end
%% Read the MAT files one by one and store the results in a cell array
% Read the concentrations and sort them in ascending order
for i=1:Nconc
    nameparts       = split(names{i},'_');
    ConcPercent(i)  = str2double(erase(nameparts{3},'w'));
end

[ConcPercent,idx]   = sort(ConcPercent,'descend');
names               = names(idx);      
VolumeData_GSB      = cell(Nconc,1);
VolumeData_ESA      = cell(Nconc,1);
StepSize            = zeros(Nconc,1);

% Now plot the stuff
for i=1:Nconc
    nameparts       = split(names{i},'_');
    PrismID(i)      = nameparts{1};
    SolutionID(i)   = nameparts{2};
    load([scriptdir filesep subdir filesep names{i}])

        %%% REDO THE NORMALISATION.
        % The volume structure dimensions are [time x peak # x (GSB/ESA)]
    switch NormType
        case 'Old'
            % Normalize the data in the "usual" way
            NormVols(:,[1 3],1) = -fitPar.Vols(:,[1 3],1)./max(abs(fitPar.Vols(:,1,1))); % GSB (Diag, Xpeak
            NormVols(:,[1 3],2) = fitPar.Vols(:,[1 3],2)./max(abs(fitPar.Vols(:,1,2)));     
            NormVols(:,[2 4],1) = -fitPar.Vols(:,[2 4],1)./max(abs(fitPar.Vols(:,2,1)));
            NormVols(:,[2 4],2) = fitPar.Vols(:,[2 4],2)./max(abs(fitPar.Vols(:,2,2)));

            NormErr(:,[1 3],1)  = fitErr.Vols(:,[1 3],1)./max(abs(fitPar.Vols(:,1,1)));
            NormErr(:,[1 3],2)  = fitErr.Vols(:,[1 3],2)./max(abs(fitPar.Vols(:,1,2)));
            NormErr(:,[2 4],1)  = fitErr.Vols(:,[2 4],1)./max(abs(fitPar.Vols(:,2,1)));
            NormErr(:,[2 4],2)  = fitErr.Vols(:,[2 4],2)./max(abs(fitPar.Vols(:,2,2)));

            NormVols(:,[3 4],:) = 10*NormVols(:,[3 4],:);   % MULTIPLY CROSS PEAKS BY 10
            NormErr(:,[3 4],:) = 10*NormErr(:,[3 4],:);     % MULTIPLY ERRORS BY 10 [ACCORDINGLY]
        case 'Sqrt'
            % Normalise by dividing by sqrt(pump*probe)
            for p=1:2 % 1=GSB, 2=ESA ---- VOLUMES
                NormVols(:,1,p) = ((-1)^p)*fitPar.Vols(:,1,p)./max(abs(fitPar.Vols(:,1,p))); % 13CO diag
                NormVols(:,2,p) = ((-1)^p)*fitPar.Vols(:,2,p)./max(abs(fitPar.Vols(:,2,p))); % 12CO diag
                NormVols(:,3,p) = ((-1)^p)*fitPar.Vols(:,3,p)./sqrt(max(abs(fitPar.Vols(:,1,p))).*max(abs(fitPar.Vols(:,2,p)))); % 13->12 uphill Xpeak
                NormVols(:,4,p) = ((-1)^p)*fitPar.Vols(:,4,p)./sqrt(max(abs(fitPar.Vols(:,1,p))).*max(abs(fitPar.Vols(:,2,p)))); % 12->13 downhill Xpeak
            end
            NormVols(:,[3 4],:) = 10*NormVols(:,[3 4],:); % MULTIPLY CROSS PEAKS BY 10
            
            for p=1:2 % 1=GSB, 2=ESA  ---- ERRORS
                NormErr(:,1,p)  = ((-1)^p)*fitErr.Vols(:,1,p)./max(abs(fitErr.Vols(:,1,p))); % 13CO diag
                NormErr(:,2,p)  = ((-1)^p)*fitErr.Vols(:,2,p)./max(abs(fitErr.Vols(:,2,p))); % 12CO diag
                NormErr(:,3,p)  = ((-1)^p)*fitErr.Vols(:,3,p)./sqrt(max(abs(fitErr.Vols(:,1,p))).*max(abs(fitErr.Vols(:,2,p)))); % 13->12 uphill Xpeak
                NormErr(:,4,p)  = ((-1)^p)*fitErr.Vols(:,4,p)./sqrt(max(abs(fitErr.Vols(:,1,p))).*max(abs(fitErr.Vols(:,2,p)))); % 12->13 downhill Xpeak
            end
            NormErr(:,[3 4],:)  = 10*NormErr(:,[3 4],:); % MULTIPLY ERRORS BY 10 [ACCORDINGLY]
        case 'Prod'
            % Normalise by dividing by (pump*probe)
            for p=1:2 % 1=GSB, 2=ESA ---- VOLUMES
                NormVols(:,1,p) = ((-1)^p)*fitPar.Vols(:,1,p)./max(abs(fitPar.Vols(:,1,p))); % 13CO diag
                NormVols(:,2,p) = ((-1)^p)*fitPar.Vols(:,2,p)./max(abs(fitPar.Vols(:,2,p))); % 12CO diag
                NormVols(:,3,p) = ((-1)^p)*fitPar.Vols(:,3,p)./(max(abs(fitPar.Vols(:,1,p))).*max(abs(fitPar.Vols(:,2,p)))); % 13->12 uphill Xpeak
                NormVols(:,4,p) = ((-1)^p)*fitPar.Vols(:,4,p)./(max(abs(fitPar.Vols(:,1,p))).*max(abs(fitPar.Vols(:,2,p)))); % 12->13 downhill Xpeak
            end
            NormVols(:,[3 4],:) = 10*NormVols(:,[3 4],:); % MULTIPLY CROSS PEAKS BY 10
            
            for p=1:2 % 1=GSB, 2=ESA  ---- ERRORS
                NormErr(:,1,p)  = ((-1)^p)*fitErr.Vols(:,1,p)./max(abs(fitErr.Vols(:,1,p))); % 13CO diag
                NormErr(:,2,p)  = ((-1)^p)*fitErr.Vols(:,2,p)./max(abs(fitErr.Vols(:,2,p))); % 12CO diag
                NormErr(:,3,p)  = ((-1)^p)*fitErr.Vols(:,3,p)./(max(abs(fitErr.Vols(:,1,p))).*max(abs(fitErr.Vols(:,2,p)))); % 13->12 uphill Xpeak
                NormErr(:,4,p)  = ((-1)^p)*fitErr.Vols(:,4,p)./(max(abs(fitErr.Vols(:,1,p))).*max(abs(fitErr.Vols(:,2,p)))); % 12->13 downhill Xpeak
            end
            NormErr(:,[3 4],:)  = 10*NormErr(:,[3 4],:); % MULTIPLY ERRORS BY 10 [ACCORDINGLY]
        case 'Square'
            % Normalise by dividing by (pump*probe)^2
            for p=1:2 % 1=GSB, 2=ESA ---- VOLUMES
                NormVols(:,1,p) = ((-1)^p)*fitPar.Vols(:,1,p)./(max(abs(fitPar.Vols(:,1,p)))).^2; % 13CO diag
                NormVols(:,2,p) = ((-1)^p)*fitPar.Vols(:,2,p)./(max(abs(fitPar.Vols(:,2,p)))).^2;    % 12CO diag
                NormVols(:,3,p) = ((-1)^p)*fitPar.Vols(:,3,p)./(max(abs(fitPar.Vols(:,1,p))).*max(abs(fitPar.Vols(:,2,p)))).^2; % 13->12 uphill Xpeak
                NormVols(:,4,p) = ((-1)^p)*fitPar.Vols(:,4,p)./(max(abs(fitPar.Vols(:,1,p))).*max(abs(fitPar.Vols(:,2,p)))).^2; % 12->13 downhill Xpeak
            end
            NormVols(:,[3 4],:) = 10*NormVols(:,[3 4],:);
            
            for p=1:2 % 1=GSB, 2=ESA  ---- ERRORS
                NormErr(:,1,p)  = ((-1)^p)*fitErr.Vols(:,1,p)./(max(abs(fitErr.Vols(:,1,p)))).^2; % 13CO diag
                NormErr(:,2,p)  = ((-1)^p)*fitErr.Vols(:,2,p)./(max(abs(fitErr.Vols(:,2,p)))).^2; % 12CO diag
                NormErr(:,3,p)  = ((-1)^p)*fitErr.Vols(:,3,p)./(max(abs(fitErr.Vols(:,1,p))).*max(abs(fitErr.Vols(:,2,p)))).^2; % 13->12 uphill Xpeak
                NormErr(:,4,p)  = ((-1)^p)*fitErr.Vols(:,4,p)./(max(abs(fitErr.Vols(:,1,p))).*max(abs(fitErr.Vols(:,2,p)))).^2; % 12->13 downhill Xpeak
            end
            NormErr(:,[3 4],:)  = 10*NormErr(:,[3 4],:);
    end
    VolumeData_GSB{i} = NormVols(:,:,1);
    VolumeData_ESA{i} = NormVols(:,:,2);
    
    writematrix(NormVols(:,:,1),[scriptdir filesep subdir filesep names{i} '_NormVols_GSB.dat']);
    writematrix(NormVols(:,:,2),[scriptdir filesep subdir filesep names{i} '_NormVols_ESA.dat']);
    
    SSR_fit(i)       = SSR;
    StepSize(i)      = output_st.stepsize;
    switch plotWhat
        case 'Diagonal GSB'
            plot(ax_UP,t2delays,NormVols(:,1,1),'^-','MarkerSize',2,'LineWidth',1.5,'Color',cmUP(i,:),'HandleVisibility','off');
            plot(ax_DW,t2delays,NormVols(:,2,1),'v-','MarkerSize',2,'LineWidth',1.5,'Color',cmDW(i,:),'HandleVisibility','off');
        case 'Xpeak GSB'
            plot(ax_UP,t2delays,NormVols(:,3,1),'^-','MarkerSize',2,'LineWidth',1.5,'Color',cmUP(i,:),'HandleVisibility','off');
            plot(ax_DW,t2delays,NormVols(:,4,1),'v-','MarkerSize',2,'LineWidth',1.5,'Color',cmDW(i,:),'HandleVisibility','off');
        case 'Diagonal ESA'
            plot(ax_UP,t2delays,NormVols(:,1,2),'^-','MarkerSize',2,'LineWidth',1.5,'Color',cmUP(i,:),'HandleVisibility','off');
            plot(ax_DW,t2delays,NormVols(:,2,2),'v-','MarkerSize',2,'LineWidth',1.5,'Color',cmDW(i,:),'HandleVisibility','off');
        case 'Xpeak ESA'
            plot(ax_UP,t2delays,NormVols(:,3,2),'^-','MarkerSize',2,'LineWidth',1.5,'Color',cmUP(i,:),'HandleVisibility','off');
            plot(ax_DW,t2delays,NormVols(:,4,2),'v-','MarkerSize',2,'LineWidth',1.5,'Color',cmDW(i,:),'HandleVisibility','off');
        case 'Xpeak GSB + C(t)'
            plot(ax_DW,t2delays,fitPar.C(:,1),'^-','MarkerSize',2,'LineWidth',1.5,'Color',cmUP(i,:),'HandleVisibility','off');
            plot(ax_UP,t2delays,NormVols(:,4,1),'v-','MarkerSize',2,'LineWidth',1.5,'Color',cmDW(i,:),'HandleVisibility','off');
            pepita(i) = fitPar.C(4,1);
        case 'Xpeak ESA + C(t)'
            plot(ax_DW,t2delays,fitPar.C(:,1),'^-','MarkerSize',2,'LineWidth',1.5,'Color',cmUP(i,:),'HandleVisibility','off');
            plot(ax_UP,t2delays,NormVols(:,4,2),'v-','MarkerSize',2,'LineWidth',1.5,'Color',cmDW(i,:),'HandleVisibility','off');
            pepita(i) = fitPar.C(4,1);
    end
    switch concType
        case '100-%'
            plot(ax_UP,NaN,NaN,'-','LineWidth',4,'Color',cmUP(i,:),'DisplayName',num2str(100-ConcPercent(i)))
            switch plotFormat
            case 'Vertical'
                plot(ax_DW,NaN,NaN,'-','LineWidth',4,'Color',cmDW(i,:),'DisplayName',num2str(100-ConcPercent(i)))
            end
        case '%'
            plot(ax_UP,NaN,NaN,'-','LineWidth',4,'Color',cmUP(i,:),'DisplayName',num2str(ConcPercent(i)))
            switch plotFormat
            case 'Vertical'
                plot(ax_DW,NaN,NaN,'-','LineWidth',4,'Color',cmDW(i,:),'DisplayName',num2str(ConcPercent(i)))
            end
    end

end

% Customize axes
xpeakUp     = 'Re(^{13}CO) \rightarrow Re(^{12}CO) \rm{\times10}';
xpeakDown   = 'Re(^{13}CO) \leftarrow Re(^{12}CO) \rm{\times10}';
titUP       = ['Uphill energy transfer' ', ' xpeakDown];
titDW       = ['Downhill energy transfer' ', ' xpeakUp];

% Set axes limits
axis(ax_UP,'tight');
axis(ax_DW,'tight');

%%% Nice formatting
% Titles
title(ax_UP,titFW,'FontSize',16);
title(ax_DW,titBW,'FontSize',16);

% Axis labels
xlabel(ax_UP,'t_2 delay (ps)','FontWeight','bold','FontSize',16);
xlabel(ax_DW,'t_2 delay (ps)','FontWeight','bold','FontSize',16);

if contains(plotWhat,'Xpeak')
    ylim(ax_UP,[-0.05 0.3]);
    ylim(ax_DW,[-0.05 0.8]);
    ylabel(ax_UP,'Normalised peak volume (\times10)','FontWeight','bold','FontSize',16);
    switch plotFormat
        case 'Vertical'
            ylabel(ax_DW,'Normalised peak volume (\times10)','FontWeight','bold','FontSize',16);
    end
else
    ylim(ax_UP,[-0.05 1]);
    ylim(ax_DW,[-0.05 1]);
    ylabel(ax_UP,'Normalised peak volume','FontWeight','bold','FontSize',16);
    switch plotFormat
        case 'Vertical'
            ylabel(ax_DW,'Normalised peak volume','FontWeight','bold','FontSize',16);
    end
end

% Axis limits
xlim(ax_UP,[0 t2delays(end)]);
xlim(ax_DW,[0 t2delays(end)]);



% Add zero line
hline_UP = yline(ax_UP,0,'HandleVisibility','off'); hline_UP.Color = [0.5 0.5 0.5];
hline_DW = yline(ax_DW,0,'HandleVisibility','off'); hline_DW.Color = [0.5 0.5 0.5];    

% Legends
switch plotFormat
    case 'Vertical'
        leg_DW = legend(ax_DW,'show');
        legend(ax_DW,'boxoff','FontWeight','bold')
        legend(ax_DW,'location','eastoutside')
        leg_UP = legend(ax_UP,'show');
        legend(ax_UP,'boxoff','FontWeight','bold')
        legend(ax_UP,'location','eastoutside')
        leg_UP.ItemTokenSize = [20,10];
        leg_DW.ItemTokenSize = [20,10];
    case 'Horizontal'
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
        ax_UP.Position(3:4) = [0.705  0.325];
        ax_DW.Position(3:4) = [0.705  0.325];
    case 'Horizontal'
        ax_UP.Position = [0.10 0.15 0.385 0.65];
        ax_DW.Position = [0.55 0.15 0.385 0.65];
end
hold(ax_UP,'off');
hold(ax_DW,'off');

if contains(plotWhat,'C(t)')
    ylim(ax_DW,[-0.05 0.8]);
    ylabel(ax_DW,'C(t_2)','FontWeight','bold','FontSize',16);
    ax_DW.Position(1) = ax_DW.Position(1)+0.025;
    title(ax_DW,'Spectral diffusion');
    title(ax_UP,'Energy transfer cross peak')
end   