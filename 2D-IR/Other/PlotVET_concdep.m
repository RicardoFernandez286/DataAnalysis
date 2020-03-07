DoFit   = 0;
DoSave  = 0;

PlotSimOnly = 0;
PlotExpOnly = 0;
PlotAllSim  = 1;
%% Get a list of MAT files
scriptdir   = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\11) 2D IR distance - na\Latest Simulations\Final Fits';
subdir      = 'SMALL';
% subdir = 'CNBz';

plotWhat    = 'Xpeak ESA'; % Xpeak or Diagonal + GSB/ESA
plotFormat  = 'Vertical'; % 'Horizontal' or 'Vertical'
concType    = '100-%'; % '100-%' or '%'
% diluent     = 'Re(^{13}C^{18}O)'; % 'Re(^{13}C^{18}O)'
diluent     = 'CNBz';
xpos        = 1.02; % 1.02 for CNBz, 0.1 horizontal
% xpos        = -0.175;
ypos        = 0.95;
FontSize    = 14;
LineWidth   = 2;
filelist    = dir([scriptdir filesep subdir]);

%% Parse the names into concentrations and make a list
names       = {filelist.name}';
names       = flipud(names(contains(names,'.mat')));
Nconc       = length(names);

% Initialise variables
ConcPercent     = zeros(Nconc,1);
PrismID         = strings(Nconc,1);
SolutionID      = strings(Nconc,1);
VolumeData_GSB  = cell(Nconc,1);
VolumeData_ESA  = cell(Nconc,1);
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
        fh.Position     = [0.3536    0.1000    0.275    0.60];
        ax_FW           = subplot(2,1,1);
        ax_BW           = subplot(2,1,2);
    case 'Horizontal'
%         fh.Position     = [0.15    0.15    0.6    0.45];
        ax_FW           = subplot(1,2,1);
        ax_BW           = subplot(1,2,2);
end

fh.Color        = [1 1 1];
fh.Units        = 'pixels';
ax_FW.FontSize  = FontSize;
ax_BW.FontSize  = FontSize;

box(ax_FW,'on');
box(ax_BW,'on');

hold(ax_FW,'on');
hold(ax_BW,'on');

IsExp= contains(names,'2D');
Nexp = sum(IsExp);
Nsimu= sum(IsExp);

if PlotAllSim == 1
    cmFW = colormap(ax_FW,(othercolor('Mrainbow',Nconc)));
    cmBW = colormap(ax_BW,(othercolor('Mrainbow',Nconc)));
else
    cmFW = colormap(ax_FW,(othercolor('Mrainbow',Nexp)));
    cmBW = colormap(ax_BW,(othercolor('Mrainbow',Nexp)));
end

% cmFW = colormap(ax_FW,flip(othercolor('Blues3',Nconc)));
% cmBW = colormap(ax_BW,flip(othercolor('Reds3',Nconc)));

%% Read the MAT files one by one and store the results in a cell array
% Read the concentrations and sort them in ascending order
for i=1:Nconc
    nameparts       = split(names{i},'_');
    if IsExp(i)
        ConcPercent(i)  = round(str2double(nameparts{3}));
    else
        ConcPercent(i)  = round(str2double(nameparts{3})*100);
    end
end

[ConcPercent,idx]   = sort(ConcPercent,'descend');
names               = names(idx);
IsExp               = IsExp(idx);

% Now plot the stuff
for i=1:Nconc
    load([scriptdir filesep subdir filesep names{i}])
    VolumeData_GSB{i} = NormVols(:,:,1);
    VolumeData_ESA{i} = NormVols(:,:,2);
    
    if PlotAllSim == 1
        plotID = i;
    else
        plotID          = find(ConcPercent(IsExp)==ConcPercent(i),1,'first');
    end
    delays{i}     = t2delays;
    xpeak_fw{i}   = NormVols(:,3,2);
    xpeak_bw{i}   = NormVols(:,4,2);
    
    if ~IsExp(i)
        if PlotExpOnly == 1
            continue
        end
        if PlotAllSim ~= 1 && ~ismember(ConcPercent(i),ConcPercent(IsExp))
            continue
        end
        decay = 1*(0.15*exp(-t2delays./2.5) + 0.85*exp(-t2delays./20));
%         decay = 1.25*(0.3*exp(-t2delays./3) + 0.7*exp(-t2delays./20)); % From JP's paper
%         decay = ones(length(t2delays),1);
        if DoFit == 0
            line_up = ':';
            line_dw = ':';
        else
            line_up = '^';
            line_dw = 'v';
        end
    else
        if PlotSimOnly == 1
            continue
        end
        decay = ones(length(t2delays),1);
        if DoFit == 0
            line_up = '^-';
            line_dw = 'v-';
        else
            line_up = '^';
            line_dw = 'v';
        end
    end

    SSR_fit(i)       = SSR;
    StepSize(i)      = output_st.stepsize;
    switch plotWhat
        case 'Diagonal GSB'
            plot(ax_FW,t2delays,NormVols(:,1,1).*decay,line_dw,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmFW(plotID,:),'HandleVisibility','off');
            plot(ax_BW,t2delays,NormVols(:,2,1).*decay,line_up,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmBW(plotID,:),'HandleVisibility','off');
        case 'Xpeak GSB'
            plot(ax_FW,t2delays,NormVols(:,3,1).*decay,line_dw,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmFW(plotID,:),'HandleVisibility','off');
            plot(ax_BW,t2delays,NormVols(:,4,1).*decay,line_up,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmBW(plotID,:),'HandleVisibility','off');
%             plot(ax_FW,t2delays,NormVols(:,3,1),'-','MarkerSize',2,'LineWidth',1.5,'Color',cmFW(i,:),'HandleVisibility','off');
%             plot(ax_BW,t2delays,NormVols(:,4,1),'-','MarkerSize',2,'LineWidth',1.5,'Color',cmBW(i,:),'HandleVisibility','off');
        case 'Diagonal ESA'
            plot(ax_FW,t2delays,NormVols(:,1,2).*decay,line_dw,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmFW(plotID,:),'HandleVisibility','off');
            plot(ax_BW,t2delays,NormVols(:,2,2).*decay,line_up,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmBW(plotID,:),'HandleVisibility','off');
        case 'Xpeak ESA'
%             plot(ax_FW,t2delays,NormVols(:,3,2)./NormVols(:,1,2).*decay,line_dw,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmFW(i,:),'HandleVisibility','off');
%             plot(ax_BW,t2delays,NormVols(:,4,2)./NormVols(:,2,2).*decay,line_up,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmBW(i,:),'HandleVisibility','off');
            plot(ax_FW,t2delays,NormVols(:,3,2).*decay,line_dw,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmFW(plotID,:),'HandleVisibility','off');
            plot(ax_BW,t2delays,NormVols(:,4,2).*decay,line_up,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmBW(plotID,:),'HandleVisibility','off');
            %             pepita(i) = NormVols(15,3,2)+NormVols(14,3,2); % WAS 3 AND 4
        case 'SpecDiff'
            plot(ax_FW,t2delays,fitPar.C(:,1),line_dw,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmFW(plotID,:),'HandleVisibility','off');
            plot(ax_BW,t2delays,fitPar.C(:,2),line_up,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmBW(plotID,:),'HandleVisibility','off');
        case 'Xpeak Diff'
            plot(ax_FW,t2delays,(NormVols(:,3,1)+NormVols(:,3,2)).*decay,line_dw,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmFW(plotID,:),'HandleVisibility','off');
            plot(ax_BW,t2delays,(NormVols(:,4,1)+NormVols(:,4,2)).*decay,line_up,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmBW(plotID,:),'HandleVisibility','off');
    end
    
    switch concType
        case '100-%'
            if ~isnan(str2double(SolutionID(i))) % Simu
%                 plot(ax_FW,NaN,NaN,'-','LineWidth',4,'Color',cmFW(plotID,:),'DisplayName',[num2str(100-ConcPercent(i)) '*'])
            else
                plot(ax_FW,NaN,NaN,'-','LineWidth',4,'Color',cmFW(plotID,:),'DisplayName',num2str(100-ConcPercent(i)))
            end
            switch plotFormat
            case 'Vertical'
                plot(ax_BW,NaN,NaN,'-','LineWidth',4,'Color',cmBW(plotID,:),'DisplayName',num2str(100-ConcPercent(i)))
            end
        case '%'
            plot(ax_FW,NaN,NaN,'-','LineWidth',4,'Color',cmFW(plotID,:),'DisplayName',num2str(ConcPercent(i)))
            switch plotFormat
            case 'Vertical'
                plot(ax_BW,NaN,NaN,'-','LineWidth',4,'Color',cmBW(plotID,:),'DisplayName',num2str(ConcPercent(i)))
            end
    end

end

% Customize axes
xpeakUp     = 'Re(^{13}CO) \rightarrow Re(^{12}CO) \rm{\times10}';
xpeakDown   = 'Re(^{13}CO) \leftarrow Re(^{12}CO) \rm{\times10}';
% titFW       = ['Downhill energy transfer' ', ' xpeakDown];
% titBW       = ['Uphill energy transfer' ', ' xpeakUp];

titFW       = 'Downhill energy transfer';
titBW       = 'Uphill energy transfer';
% Set axes limits
axis(ax_FW,'tight');
axis(ax_BW,'tight');

%%% Nice formatting
% Titles
title(ax_FW,titFW,'FontSize',FontSize);
title(ax_BW,titBW,'FontSize',FontSize);

% Axis labels
xlabel(ax_FW,'t_2 delay (ps)','FontWeight','bold','FontSize',FontSize);
xlabel(ax_BW,'t_2 delay (ps)','FontWeight','bold','FontSize',FontSize);

if contains(plotWhat,'Xpeak')
    ylim(ax_FW,[-0.05 0.8]);
    ylim(ax_BW,[-0.05 0.4]);
    ylabel(ax_FW,'Norm. peak volume (\times10)','FontWeight','bold','FontSize',FontSize);
    switch plotFormat
        case 'Vertical'
            ylabel(ax_BW,'Norm. peak volume (\times10)','FontWeight','bold','FontSize',FontSize);
    end
elseif contains(plotWhat,'SpecDiff')
    ylim(ax_FW,[-0.05 1]);
    ylim(ax_BW,[-0.05 1]);
    ylabel(ax_FW,'C (t_{2})','FontWeight','bold','FontSize',FontSize);
    switch plotFormat
        case 'Vertical'
            ylabel(ax_BW,'C (t_{2})','FontWeight','bold','FontSize',FontSize);
    end
else
    ylim(ax_FW,[-0.05 1]);
    ylim(ax_BW,[-0.05 1]);
    ylabel(ax_FW,'Normalised peak volume','FontWeight','bold','FontSize',FontSize);
    switch plotFormat
        case 'Vertical'
            ylabel(ax_BW,'Normalised peak volume','FontWeight','bold','FontSize',FontSize);
    end
end


% Axis limits
% xlim(ax_FW,[0 t2delays(end)]);
% xlim(ax_BW,[0 t2delays(end)]);

xlim(ax_FW,[0 60]);
xlim(ax_BW,[0 60]);

% ylim(ax_FW,[-0.05 1]);
% ylim(ax_BW,[-0.05 1]);
axis(ax_FW,'tight')
axis(ax_BW,'tight')

% Add zero line
hline_FW = yline(ax_FW,0,'HandleVisibility','off'); hline_FW.Color = [0.5 0.5 0.5];
hline_BW = yline(ax_BW,0,'HandleVisibility','off'); hline_BW.Color = [0.5 0.5 0.5];

% Legends
switch plotFormat
    case 'Vertical'
        text(ax_FW,xpos,ypos,['\bf{% ' diluent '}'],'FontSize',12,'FontWeight','bold','Units','normalized')
        plot(ax_FW,NaN,NaN,'.','Color','w','DisplayName','')
%         text(ax_BW,xpos,ypos,['\bf{% ' diluent '}'],'FontSize',12,'FontWeight','bold','Units','normalized')
        plot(ax_BW,NaN,NaN,'.','Color','w','DisplayName','')
        leg_BW = legend(ax_BW,'show');
        legend(ax_BW,'boxoff','FontWeight','bold')
        legend(ax_BW,'location','eastoutside')
        leg_BW.ItemTokenSize = [20,10];
        leg_FW = legend(ax_FW,'show');
        legend(ax_FW,'boxoff','FontWeight','bold')
        legend(ax_FW,'location','bestoutside')
        leg_FW.ItemTokenSize = [20,10];
        delete(leg_BW);
    case 'Horizontal'
        text(ax_FW,xpos,1.2,['\bf{% ' diluent '}'],'FontSize',12,'FontWeight','bold','Units','normalized')
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


ax_FW.XTick = [0:10:max(t2delays)];
ax_BW.XTick = [0:10:max(t2delays)];
ax_FW.TickLength = 1.6*ax_FW.TickLength;
ax_BW.TickLength = 1.6*ax_BW.TickLength;
leg_FW.Position(2) = leg_FW.Position(2)-0.04;

%%%% INTEGRATE THE TIME-DEPENTENT KINETICS
integral_bw = zeros(Nconc,1);
integral_fw = zeros(Nconc,1);

for i=1:Nconc
    integral_bw(i) = trapz(delays{i},xpeak_bw{i});
    integral_fw(i) = trapz(delays{i},xpeak_fw{i});
end

if DoSave == 1
    for i=1:Nconc
        ESA_down(:,i)    = VolumeData_ESA{i}(:,3);
        ESA_up(:,i)      = VolumeData_ESA{i}(:,4);
    end
    dlmwrite([scriptdir filesep subdir 'ESA_DOWN.dat'],[[0 ConcPercent'];[t2delays ESA_down.*decay]]);
    dlmwrite([scriptdir filesep subdir 'ESA_UP.dat'],[[0 ConcPercent'];[t2delays ESA_up.*decay]]);
end

if DoFit == 0
 return
end

%% Do kinetic model fit
for i=1:Nconc
    FW_data = [VolumeData_ESA{i}(:,1) VolumeData_ESA{i}(:,3)/10];
    BW_data = [VolumeData_ESA{i}(:,2) VolumeData_ESA{i}(:,4)/10];
    time = delays{i};

    fitPar0 = [3  20   70   90  0.3];
    LB      = [2  0    0    0   0  ];
    UB      = [5  50   500  500 1  ];

    ShowOutput = 0;
    DoFit = 1;
    tic
    fitpar(:,i) = Fit_EnT(time,FW_data,BW_data,fitPar0,LB,UB,ShowOutput,DoFit,fh);
    disp(['Fitting done for ' num2str(i) ' of ' num2str(Nconc) ' (' num2str(toc) 's)']);
end

fitpar = fitpar';


% FIGURE SIZE
fh.Units = 'pixels';
% fh.Position = [275 718 1285 620];