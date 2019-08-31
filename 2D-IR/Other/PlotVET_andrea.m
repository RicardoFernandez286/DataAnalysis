%% Get a list of MAT files
scriptdir   = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\9) 2D IR distance - na\Data\New fits';
subdir      = 'Andrea';

plotWhat    = 'Xpeak GSB + C(t)'; 
% plotWhat    = 'Xpeak ESA'; % Xpeak or Diagonal + GSB/ESA
plotFormat  = 'Horizontal'; % 'Horizontal' or 'Vertical'
concType    = '%'; % '100-%' or '%'
diluent     = 'D_{2}O'; % 'Re(^{13}C^{18}O)'

% diluent     = 'Re(^{13}C^{18}O)';
xpos        = 0; % 1.02 for CNBz, 0.95 for Re18, 0.85 for D2O
                % 0.25 for horizontal

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
        ax_FW           = subplot(2,1,1);
        ax_BW           = subplot(2,1,2);
    case 'Horizontal'
        fh.Position     = [0.15    0.15    0.6    0.45];
        ax_FW           = subplot(1,2,1);
        ax_BW           = subplot(1,2,2);
end

fh.Color        = [1 1 1];
fh.Units        = 'pixels';
ax_FW.FontSize  = 16;
ax_BW.FontSize  = 16;

box(ax_FW,'on');
box(ax_BW,'on');

hold(ax_FW,'on');
hold(ax_BW,'on');

cmFW = colormap(ax_FW,othercolor('Mrainbow',Nconc));
cmBW = colormap(ax_BW,othercolor('Mrainbow',Nconc));

% cmFW = colormap(ax_FW,flip(othercolor('Blues3',Nconc)));
% cmBW = colormap(ax_BW,flip(othercolor('Reds3',Nconc)));

switch plotFormat
    case 'Vertical'
        text(ax_BW,xpos,0.94,['\bf{% ' diluent '}'],'FontSize',12,'FontWeight','bold','Units','normalized')
        plot(ax_BW,NaN,NaN,'.','Color','w','DisplayName','')
        text(ax_FW,xpos,0.94,['\bf{% ' diluent '}'],'FontSize',12,'FontWeight','bold','Units','normalized')
        plot(ax_FW,NaN,NaN,'.','Color','w','DisplayName','')
    case 'Horizontal'
        text(ax_FW,xpos,1.19,['\bf{% ' diluent '}'],'FontSize',12,'FontWeight','bold','Units','normalized')
end
%% Read the MAT files one by one and store the results in a cell array
% Read the concentrations and sort them in ascending order
for i=1:Nconc
    nameparts       = split(names{i},'_');
    ConcPercent(i)  = str2double(erase(nameparts{4},'w'));
end

[ConcPercent,idx]   = sort(ConcPercent,'descend');
names               = names(idx);      

% Now plot the stuff
for i=1:Nconc
    nameparts       = split(names{i},'_');
    PrismID(i)      = nameparts{1};
    SolutionID(i)   = nameparts{2};
    load([scriptdir filesep subdir filesep names{i}])
    VoumeData_GSB{i} = NormVols(:,:,1);
    VoumeData_ESA{i} = NormVols(:,:,2);
    SSR_fit(i)       = SSR;
    StepSize(i)      = output_st.stepsize;
    switch plotWhat
        case 'Diagonal GSB'
            plot(ax_FW,t2delays,NormVols(:,1,1),'^-','MarkerSize',2,'LineWidth',1.5,'Color',cmFW(i,:),'HandleVisibility','off');
            plot(ax_BW,t2delays,NormVols(:,2,1),'v-','MarkerSize',2,'LineWidth',1.5,'Color',cmBW(i,:),'HandleVisibility','off');
        case 'Xpeak GSB'
            plot(ax_FW,t2delays,NormVols(:,3,1),'^-','MarkerSize',2,'LineWidth',1.5,'Color',cmFW(i,:),'HandleVisibility','off');
            plot(ax_BW,t2delays,NormVols(:,4,1),'v-','MarkerSize',2,'LineWidth',1.5,'Color',cmBW(i,:),'HandleVisibility','off');
        case 'Diagonal ESA'
            plot(ax_FW,t2delays,NormVols(:,1,2),'^-','MarkerSize',2,'LineWidth',1.5,'Color',cmFW(i,:),'HandleVisibility','off');
            plot(ax_BW,t2delays,NormVols(:,2,2),'v-','MarkerSize',2,'LineWidth',1.5,'Color',cmBW(i,:),'HandleVisibility','off');
        case 'Xpeak ESA'
            plot(ax_FW,t2delays,NormVols(:,3,2),'^-','MarkerSize',2,'LineWidth',1.5,'Color',cmFW(i,:),'HandleVisibility','off');
            plot(ax_BW,t2delays,NormVols(:,4,2),'v-','MarkerSize',2,'LineWidth',1.5,'Color',cmBW(i,:),'HandleVisibility','off');
        case 'Xpeak GSB + C(t)'
            plot(ax_BW,t2delays,fitPar.C(:,1),'^-','MarkerSize',2,'LineWidth',1.5,'Color',cmFW(i,:),'HandleVisibility','off');
            plot(ax_FW,t2delays,NormVols(:,4,1),'v-','MarkerSize',2,'LineWidth',1.5,'Color',cmBW(i,:),'HandleVisibility','off');
            pepita(i) = fitPar.C(4,1);
        case 'Xpeak ESA + C(t)'
            plot(ax_BW,t2delays,fitPar.C(:,1),'^-','MarkerSize',2,'LineWidth',1.5,'Color',cmFW(i,:),'HandleVisibility','off');
            plot(ax_FW,t2delays,NormVols(:,4,2),'v-','MarkerSize',2,'LineWidth',1.5,'Color',cmBW(i,:),'HandleVisibility','off');
            pepita(i) = fitPar.C(4,1);
    end
    switch concType
        case '100-%'
            plot(ax_FW,NaN,NaN,'-','LineWidth',4,'Color',cmFW(i,:),'DisplayName',num2str(100-ConcPercent(i)))
            switch plotFormat
            case 'Vertical'
                plot(ax_BW,NaN,NaN,'-','LineWidth',4,'Color',cmBW(i,:),'DisplayName',num2str(100-ConcPercent(i)))
            end
        case '%'
            plot(ax_FW,NaN,NaN,'-','LineWidth',4,'Color',cmFW(i,:),'DisplayName',num2str(ConcPercent(i)))
            switch plotFormat
            case 'Vertical'
                plot(ax_BW,NaN,NaN,'-','LineWidth',4,'Color',cmBW(i,:),'DisplayName',num2str(ConcPercent(i)))
            end
    end

end

% Customize axes
xpeakUp     = 'Re(^{13}CO) \rightarrow Re(^{12}CO) \rm{\times10}';
xpeakDown   = 'Re(^{13}CO) \leftarrow Re(^{12}CO) \rm{\times10}';
titFW       = 'Uphill energy transfer';
titBW       = 'Downhill energy transfer';

% Set axes limits
axis(ax_FW,'tight');
axis(ax_BW,'tight');

%%% Nice formatting
% Titles
title(ax_FW,titFW,'FontSize',16);
title(ax_BW,titBW,'FontSize',16);

% Axis labels
xlabel(ax_FW,'t_2 delay (ps)','FontWeight','bold','FontSize',16);
xlabel(ax_BW,'t_2 delay (ps)','FontWeight','bold','FontSize',16);

if contains(plotWhat,'Xpeak')
    ylim(ax_FW,[-0.05 0.3]);
    ylim(ax_BW,[-0.05 0.8]);
    ylabel(ax_FW,'Normalised peak volume (\times10)','FontWeight','bold','FontSize',16);
    switch plotFormat
        case 'Vertical'
            ylabel(ax_BW,'Normalised peak volume (\times10)','FontWeight','bold','FontSize',16);
    end
else
    ylim(ax_FW,[-0.05 1]);
    ylim(ax_BW,[-0.05 1]);
    ylabel(ax_FW,'Normalised peak volume','FontWeight','bold','FontSize',16);
    switch plotFormat
        case 'Vertical'
            ylabel(ax_BW,'Normalised peak volume','FontWeight','bold','FontSize',16);
    end
end

% Axis limits
xlim(ax_FW,[0 t2delays(end)]);
xlim(ax_BW,[0 t2delays(end)]);



% Add zero line
hline_FW = yline(ax_FW,0,'HandleVisibility','off'); hline_FW.Color = [0.5 0.5 0.5];
hline_BW = yline(ax_BW,0,'HandleVisibility','off'); hline_BW.Color = [0.5 0.5 0.5];    

% Legends
switch plotFormat
    case 'Vertical'
        leg_BW = legend(ax_BW,'show');
        legend(ax_BW,'boxoff','FontWeight','bold')
        legend(ax_BW,'location','eastoutside')
        leg_FW = legend(ax_FW,'show');
        legend(ax_FW,'boxoff','FontWeight','bold')
        legend(ax_FW,'location','eastoutside')
        leg_FW.ItemTokenSize = [20,10];
        leg_BW.ItemTokenSize = [20,10];
    case 'Horizontal'
        leg_FW = legend(ax_FW,'show');
        legend(ax_FW,'boxoff')
        legend(ax_FW,'FontWeight','bold')
        legend(ax_FW,'orientation','horizontal')
        leg_FW.Position = [0.3 0.9 0.5 0.06];
        leg_FW.ItemTokenSize = [20,10];
end

% Make figure resizable
fh.Units        = 'normalized';

switch plotFormat
    case 'Vertical'
        ax_FW.Position(3:4) = [0.705  0.325];
        ax_BW.Position(3:4) = [0.705  0.325];
    case 'Horizontal'
        ax_FW.Position = [0.10 0.15 0.385 0.65];
        ax_BW.Position = [0.55 0.15 0.385 0.65];
end
hold(ax_FW,'off');
hold(ax_BW,'off');

if contains(plotWhat,'C(t)')
    ylim(ax_BW,[-0.05 0.8]);
    ylabel(ax_BW,'C(t_2)','FontWeight','bold','FontSize',16);
    ax_BW.Position(1) = ax_BW.Position(1)+0.025;
    title(ax_BW,'Spectral diffusion');
    title(ax_FW,'Energy transfer cross peak')
end   