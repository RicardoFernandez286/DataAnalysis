DoFit   = 0;
DoSave  = 0;

NormType='Sqrt'; % 'Old' 'Sqrt' 'Prod' 'Squared'

% PlotWhat='Pick';
% PlotWhat='E';
PlotWhat='Exp';
% PlotWhat='Sim';
% PlotWhat='Mixed';

PlotSimOnly = 0;
PlotExpOnly = 0;
PlotAllSim  = 1;
%% Get a list of MAT files
switch PlotWhat
    case 'E'
        scriptdir = 'E:\';
        subdir = 'FitResults';
    case 'Sim'
%         scriptdir   = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\Submitted\11) VET distance - JPC\OLD\Latest Simulations\Final simulations';       
%         subdir      = '\SmallDil_BiggerBox_NOClustering\FitResults';
%         subdir      = '\SmallDil_BiggerBox_Clustering\FitResults';
%         subdir      = 'CNBz';
        scriptdir   = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\Submitted\18) VET distance - JPC\OLD\Data\Final Fits\Simulations';
        subdir      = 'Big molecule';
%         subdir      = 'Small molecule - NO clustering';
%         subdir      = 'Small molecule - WITH clustering';
%         scriptdir = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\Submitted\18) VET distance - JPC\OLD\Latest Simulations\Final simulations\Development version - cone\For Polaron\Cone_dil\FitResults';
%         subdir = [];
    case 'Exp'
        scriptdir   = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\Submitted\18) VET distance - JPC\OLD\Data\Final Fits\Experiment';
        subdir      = 'CNBz - All';
%         subdir      = 'Re18 - All';
%         subdir = [];
    case 'Pick'
        scriptdir = uigetdir();
        subdir = [];
    case 'Mixed'
        scriptdir   = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\Submitted\18) VET distance - JPC\OLD\Data\Final Fits\';
        subdir      = 'Mixed\Big';
end
plotWhat    = 'Xpeak ESA'; % Xpeak or Diagonal + GSB/ESA
plotFormat  = 'Vertical'; % 'Horizontal' or 'Vertical'
concType    = '%'; % '100-%' or '%'
% diluent     = 'Re(^{13}C^{18}O)'; % 'Re(^{13}C^{18}O)'
diluent     = 'Re';

xpos        = 0.83; % 0.8 for Re18; 0.83 for CNBz

% xpos        = 1.04; % 1.04 for Re18; 1.06 for CNBz; 0.1 horizontal
% xpos        = -0.175;
ypos        = 0.35;

FontSize    = 14;
LineWidth   = 2;
filelist    = dir([scriptdir filesep subdir]);

%% Parse the names into concentrations and make a list
names       = {filelist.name}';
names       = flipud(names(contains(names,'.mat')));
% names       = names(1:end-1);
Nconc       = length(names);
dilupos     = 3;

% Initialise variables
ConcPercent     = zeros(Nconc,1);
PrismID         = strings(Nconc,1);
SolutionID      = strings(Nconc,1);
VolumeData_GSB  = cell(Nconc,1);
VolumeData_ESA  = cell(Nconc,1);
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
        fh.Position     = [0.3536    0.1000    0.275    0.60];
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

IsExp= contains(names,'2D');
Nexp = sum(IsExp);
Nsimu= sum(IsExp);

if PlotAllSim == 1
    cmFW = colormap(ax_UP,(othercolor('Mrainbow',Nconc)));
    cmBW = colormap(ax_DW,(othercolor('Mrainbow',Nconc)));
else
    cmFW = colormap(ax_UP,(othercolor('Mrainbow',Nexp)));
    cmBW = colormap(ax_DW,(othercolor('Mrainbow',Nexp)));
end

% cmFW = colormap(ax_UP,flip(othercolor('Blues3',Nconc)));
% cmBW = colormap(ax_DW,flip(othercolor('Reds3',Nconc)));

%% Read the MAT files one by one and store the results in a cell array
% Read the concentrations and sort them in ascending order
for i=1:Nconc
    nameparts       = split(names{i},'_');
    if IsExp(i)
        ConcPercent(i)  = round(str2double(nameparts{dilupos}));
    else
        if ~isnan(str2double(nameparts{dilupos}))
            ConcPercent(i)  = round(str2double(nameparts{dilupos})*100);
        else
            ConcPercent(i)  = round(str2double(nameparts{dilupos-1})*100);
        end
    end
end

[ConcPercent,idx]   = sort(ConcPercent,'descend');
names               = names(idx);
IsExp               = IsExp(idx);

% Now plot the stuff
for i=1:Nconc
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
        case 'Sqrt'
            % Normalise by dividing by sqrt(pump*probe)
            for p=1:2 % 1=GSB, 2=ESA ---- VOLUMES
                NormVols(:,1,p) = ((-1)^p)*fitPar.Vols(:,1,p)./max(abs(fitPar.Vols(:,1,p))); % 13CO diag
                NormVols(:,2,p) = ((-1)^p)*fitPar.Vols(:,2,p)./max(abs(fitPar.Vols(:,2,p))); % 12CO diag
                NormVols(:,3,p) = ((-1)^p)*fitPar.Vols(:,3,p)./sqrt(max(abs(fitPar.Vols(:,1,p))).*max(abs(fitPar.Vols(:,2,p)))); % 13->12 uphill Xpeak
                NormVols(:,4,p) = ((-1)^p)*fitPar.Vols(:,4,p)./sqrt(max(abs(fitPar.Vols(:,1,p))).*max(abs(fitPar.Vols(:,2,p)))); % 12->13 downhill Xpeak
            end
%             NormVols(:,[3 4],:) = 10*NormVols(:,[3 4],:);
           
            for p=1:2 % 1=GSB, 2=ESA  ---- ERRORS
                NormErr(:,1,p)  = ((-1)^p)*fitErr.Vols(:,1,p)./max(abs(fitErr.Vols(:,1,p))); % 13CO diag
                NormErr(:,2,p)  = ((-1)^p)*fitErr.Vols(:,2,p)./max(abs(fitErr.Vols(:,2,p))); % 12CO diag
                NormErr(:,3,p)  = ((-1)^p)*fitErr.Vols(:,3,p)./sqrt(max(abs(fitErr.Vols(:,1,p))).*max(abs(fitErr.Vols(:,2,p)))); % 13->12 uphill Xpeak
                NormErr(:,4,p)  = ((-1)^p)*fitErr.Vols(:,4,p)./sqrt(max(abs(fitErr.Vols(:,1,p))).*max(abs(fitErr.Vols(:,2,p)))); % 12->13 downhill Xpeak
            end
%             NormErr(:,[3 4],:)  = 10*NormErr(:,[3 4],:);                
        case 'Prod'
            % Normalise by dividing by (pump*probe)
            for p=1:2 % 1=GSB, 2=ESA ---- VOLUMES
                NormVols(:,1,p) = ((-1)^p)*fitPar.Vols(:,1,p)./max(abs(fitPar.Vols(:,1,p))); % 13CO diag
                NormVols(:,2,p) = ((-1)^p)*fitPar.Vols(:,2,p)./max(abs(fitPar.Vols(:,2,p))); % 12CO diag
                NormVols(:,3,p) = ((-1)^p)*fitPar.Vols(:,3,p)./(max(abs(fitPar.Vols(:,1,p))).*max(abs(fitPar.Vols(:,2,p)))); % 13->12 uphill Xpeak
                NormVols(:,4,p) = ((-1)^p)*fitPar.Vols(:,4,p)./(max(abs(fitPar.Vols(:,1,p))).*max(abs(fitPar.Vols(:,2,p)))); % 12->13 downhill Xpeak
            end
            
            for p=1:2 % 1=GSB, 2=ESA  ---- ERRORS
                NormErr(:,1,p)  = ((-1)^p)*fitErr.Vols(:,1,p)./max(abs(fitErr.Vols(:,1,p))); % 13CO diag
                NormErr(:,2,p)  = ((-1)^p)*fitErr.Vols(:,2,p)./max(abs(fitErr.Vols(:,2,p))); % 12CO diag
                NormErr(:,3,p)  = ((-1)^p)*fitErr.Vols(:,3,p)./(max(abs(fitErr.Vols(:,1,p))).*max(abs(fitErr.Vols(:,2,p)))); % 13->12 uphill Xpeak
                NormErr(:,4,p)  = ((-1)^p)*fitErr.Vols(:,4,p)./(max(abs(fitErr.Vols(:,1,p))).*max(abs(fitErr.Vols(:,2,p)))); % 12->13 downhill Xpeak
            end
        case 'Square'
            % Normalise by dividing by (pump*probe)^2
            for p=1:2 % 1=GSB, 2=ESA ---- VOLUMES
                NormVols(:,1,p) = ((-1)^p)*fitPar.Vols(:,1,p)./(max(abs(fitPar.Vols(:,1,p)))).^2; % 13CO diag
                NormVols(:,2,p) = ((-1)^p)*fitPar.Vols(:,2,p)./(max(abs(fitPar.Vols(:,2,p)))).^2;    % 12CO diag
                NormVols(:,3,p) = ((-1)^p)*fitPar.Vols(:,3,p)./(max(abs(fitPar.Vols(:,1,p))).*max(abs(fitPar.Vols(:,2,p)))).^2; % 13->12 uphill Xpeak
                NormVols(:,4,p) = ((-1)^p)*fitPar.Vols(:,4,p)./(max(abs(fitPar.Vols(:,1,p))).*max(abs(fitPar.Vols(:,2,p)))).^2; % 12->13 downhill Xpeak
            end
            
            for p=1:2 % 1=GSB, 2=ESA  ---- ERRORS
                NormErr(:,1,p)  = ((-1)^p)*fitErr.Vols(:,1,p)./(max(abs(fitErr.Vols(:,1,p)))).^2; % 13CO diag
                NormErr(:,2,p)  = ((-1)^p)*fitErr.Vols(:,2,p)./(max(abs(fitErr.Vols(:,2,p)))).^2; % 12CO diag
                NormErr(:,3,p)  = ((-1)^p)*fitErr.Vols(:,3,p)./(max(abs(fitErr.Vols(:,1,p))).*max(abs(fitErr.Vols(:,2,p)))).^2; % 13->12 uphill Xpeak
                NormErr(:,4,p)  = ((-1)^p)*fitErr.Vols(:,4,p)./(max(abs(fitErr.Vols(:,1,p))).*max(abs(fitErr.Vols(:,2,p)))).^2; % 12->13 downhill Xpeak
            end
    end
    NormVols    = 100*NormVols; % Now in % population
    NormErr     = 100*NormErr; % Now in % population
%     
%     NormVols(:,3,2) = NormVols(:,3,2)*45/mean(NormVols((end-2):end,3,2));
%     NormVols(:,4,2) = NormVols(:,4,2)*48/mean(NormVols((end-2):end,4,2));
%     
    VolumeData_GSB{i} = NormVols(:,:,1);
    VolumeData_ESA{i} = NormVols(:,:,2);
    
    X0Pos(:,i) = fitPar.X0;
    AMPRATIO(i,:) = squeeze(fitPar.Amps(1,1:2,2).*sqrt(fitPar.Sx(1:2))'.*sqrt(fitPar.Sy(1:2))');
    
    AmpQuot = AMPRATIO(:,2)./AMPRATIO(:,1);
    
    if ConcPercent(i) == 20
        NormVols(:,4,:) = NormVols(:,4,:);
    end

 %%   
    if PlotAllSim == 1
        plotID = i;
    else
        plotID          = find(ConcPercent(IsExp)==ConcPercent(i),1,'first');
    end
   
    if ~IsExp(i)
        if PlotExpOnly == 1
            continue
        end
        if PlotAllSim ~= 1 && ~ismember(ConcPercent(i),ConcPercent(IsExp))
            continue
        end
        a = 1;
%         decay = 1*((1-a)*exp(-t2delays./2.5) + a*exp(-t2delays./20));
        decay = exp(-t2delays./20);
%         decay = ones(length(t2delays),1);
        if DoFit == 0
            if contains(names{i},'Clus')
                line_up = '-';
                line_dw = '-';
            else
                line_up = '-';
                line_dw = '-';
            end
        else
            line_up = '-';
            line_dw = '-';
        end
    else
        if PlotSimOnly == 1
            continue
        end
        decay = ones(length(t2delays),1);
        if DoFit == 0
            line_up = '-^';
            line_dw = '-v';
        else
            line_up = '-^';
            line_dw = '-v';
        end
    end
    delays{i}     = t2delays;
    xpeak_UP{i}   = NormVols(:,3,2).*decay;
    xpeak_DW{i}   = NormVols(:,4,2).*decay;
    
    SSR_fit(i)       = SSR;
    StepSize(i)      = output_st.stepsize;
    switch plotWhat
        case 'Diagonal GSB'
            plot(ax_UP,t2delays,NormVols(:,1,1).*decay,line_up,'MarkerSize',3,'LineWidth',LineWidth,'Color',cmFW(plotID,:),'HandleVisibility','off');
            plot(ax_DW,t2delays,NormVols(:,2,1).*decay,line_dw,'MarkerSize',3,'LineWidth',LineWidth,'Color',cmBW(plotID,:),'HandleVisibility','off');
        case 'Xpeak GSB'
            plot(ax_UP,t2delays,NormVols(:,3,1).*decay,line_up,'MarkerSize',3,'LineWidth',LineWidth,'Color',cmFW(plotID,:),'HandleVisibility','off');
            plot(ax_DW,t2delays,NormVols(:,4,1).*decay,line_dw,'MarkerSize',3,'LineWidth',LineWidth,'Color',cmBW(plotID,:),'HandleVisibility','off');
%             plot(ax_UP,t2delays,NormVols(:,3,1),'-','MarkerSize',2,'LineWidth',1.5,'Color',cmFW(i,:),'HandleVisibility','off');
%             plot(ax_DW,t2delays,NormVols(:,4,1),'-','MarkerSize',2,'LineWidth',1.5,'Color',cmBW(i,:),'HandleVisibility','off');
        case 'Diagonal ESA'
            plot(ax_UP,t2delays,NormVols(:,1,2).*decay,line_up,'MarkerSize',3,'LineWidth',LineWidth,'Color',cmFW(plotID,:),'HandleVisibility','off');
            plot(ax_DW,t2delays,NormVols(:,2,2).*decay,line_dw,'MarkerSize',3,'LineWidth',LineWidth,'Color',cmBW(plotID,:),'HandleVisibility','off');
        case 'Xpeak ESA'
%             plot(ax_UP,t2delays,NormVols(:,3,2)./NormVols(:,1,2).*decay,line_up,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmFW(i,:),'HandleVisibility','off');
%             plot(ax_DW,t2delays,NormVols(:,4,2)./NormVols(:,2,2).*decay,line_dw,'MarkerSize',2,'LineWidth',LineWidth,'Color',cmBW(i,:),'HandleVisibility','off');
            plot(ax_UP,t2delays,NormVols(:,3,2).*decay,line_up,'MarkerSize',3,'LineWidth',LineWidth,'Color',cmFW(plotID,:),'HandleVisibility','off');
            plot(ax_DW,t2delays,NormVols(:,4,2).*decay,line_dw,'MarkerSize',3,'LineWidth',LineWidth,'Color',cmBW(plotID,:),'HandleVisibility','off');
            %             pepita(i) = NormVols(15,3,2)+NormVols(14,3,2); % WAS 3 AND 4
        case 'SpecDiff'
            plot(ax_UP,t2delays,fitPar.C(:,1),line_up,'MarkerSize',3,'LineWidth',LineWidth,'Color',cmFW(plotID,:),'HandleVisibility','off');
            plot(ax_DW,t2delays,fitPar.C(:,2),line_dw,'MarkerSize',3,'LineWidth',LineWidth,'Color',cmBW(plotID,:),'HandleVisibility','off');
        case 'Xpeak Diff'
            plot(ax_UP,t2delays,(NormVols(:,3,1)+NormVols(:,3,2)).*decay,line_up,'MarkerSize',3,'LineWidth',LineWidth,'Color',cmFW(plotID,:),'HandleVisibility','off');
            plot(ax_DW,t2delays,(NormVols(:,4,1)+NormVols(:,4,2)).*decay,line_dw,'MarkerSize',3,'LineWidth',LineWidth,'Color',cmBW(plotID,:),'HandleVisibility','off');
    end
    
    switch concType
        case '100-%'
            if ~isnan(str2double(SolutionID(i))) % Simu
%                 plot(ax_UP,NaN,NaN,'-','LineWidth',4,'Color',cmFW(plotID,:),'DisplayName',[num2str(100-ConcPercent(i)) '*'])
            else
                plot(ax_UP,NaN,NaN,'-','LineWidth',4,'Color',cmFW(plotID,:),'DisplayName',num2str(100-ConcPercent(i)))
            end
            switch plotFormat
            case 'Vertical'
                plot(ax_DW,NaN,NaN,'-','LineWidth',4,'Color',cmBW(plotID,:),'DisplayName',num2str(100-ConcPercent(i)))
            end
        case '%'
            plot(ax_UP,NaN,NaN,'-','LineWidth',4,'Color',cmFW(plotID,:),'DisplayName',num2str(ConcPercent(i)))
            switch plotFormat
            case 'Vertical'
                plot(ax_DW,NaN,NaN,'-','LineWidth',4,'Color',cmBW(plotID,:),'DisplayName',num2str(ConcPercent(i)))
            end
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
%     ylabel(ax_UP,'Norm. peak volume (\times10)','FontWeight','bold','FontSize',FontSize);
    ylabel(ax_UP,'Population (%)','FontWeight','bold','FontSize',FontSize);
    switch plotFormat
        case 'Vertical'
            ylabel(ax_DW,'Norm. peak volume (\times10)','FontWeight','bold','FontSize',FontSize);
            ylabel(ax_DW,'Population (%)','FontWeight','bold','FontSize',FontSize);
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
%     ylabel(ax_UP,'Normalised peak volume','FontWeight','bold','FontSize',FontSize);
    ylabel(ax_UP,'Population (%)','FontWeight','bold','FontSize',FontSize);
    switch plotFormat
        case 'Vertical'
%             ylabel(ax_DW,'Normalised peak volume','FontWeight','bold','FontSize',FontSize);
            ylabel(ax_DW,'Population (%)','FontWeight','bold','FontSize',FontSize);
    end
end


% Axis limits
% xlim(ax_UP,[0 t2delays(end)]);
% xlim(ax_DW,[0 t2delays(end)]);

% ylim(ax_UP,[-0.05 1]);
% ylim(ax_DW,[-0.05 1]);
axis(ax_UP,'tight');
axis(ax_DW,'tight');

xlim(ax_UP,[0 60]);
xlim(ax_DW,[0 60]);

% xlim(ax_UP,[0 10000]);
% xlim(ax_DW,[0 10000]);

% xlim(ax_UP,[0 max(t2delays)]);
% xlim(ax_DW,[0 max(t2delays)]);

ylim(ax_UP,[0 6]);
ylim(ax_DW,[0 6]);

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
        title(ax_UP,'(A) Uphill transfer (Peak 1)','Units','normalized','position',[0.3 1.05]);
        title(ax_DW,'(B) Downhill transfer (Peak 2)','Units','normalized','position',[0.33 1.05]);
%         title(ax_UP,'(A) ^{13}CO diagonal peak','Units','normalized','position',[0.25 1.05]);
%         title(ax_DW,'(B) ^{12}CO diagonal peak','Units','normalized','position',[0.25 1.05]);
        leg_UP.Position(1:2) = [0.85  (0.5-leg_UP.Position(4)/2)];
        annotation(fh,'textbox',[xpos -0.05+leg_UP.Position(4)+leg_UP.Position(2) 0.5 0.05],'String',['\bf{% ' diluent '}'],'FontSize',12,'FontWeight','bold','Units','normalized','EdgeColor','none');
    case 'Horizontal'
        ax_UP.Position = [0.10 0.15 0.385 0.65];
        ax_DW.Position = [0.55 0.15 0.385 0.65];
end
hold(ax_UP,'off');
hold(ax_DW,'off');



% ax_UP.XTick = 0:200:max(t2delays);
% ax_DW.XTick = 0:200:max(t2delays);

ax_UP.XScale = 'lin';
ax_DW.XScale = 'lin';

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

%% INTEGRATE THE TIME-DEPENTENT KINETICS
integral_DW = zeros(Nconc,1);
integral_UP = zeros(Nconc,1);

for i=1:Nconc
    t2delays_all(:,i) = delays{i};
    integral_DW(i) = trapz(delays{i},xpeak_DW{i});
    integral_UP(i) = trapz(delays{i},xpeak_UP{i});
end
%%
if DoSave == 1
    for i=1:Nconc
        ESA_down(:,i)    = VolumeData_ESA{i}(:,3);
        ESA_up(:,i)      = VolumeData_ESA{i}(:,4);
    end
    dlmwrite([scriptdir filesep subdir 'ESA_DOWN.dat'],[[0 ConcPercent'];[t2delays ESA_down.*decay]]);
    dlmwrite([scriptdir filesep subdir 'ESA_UP.dat'],[[0 ConcPercent'];[t2delays ESA_up.*decay]]);
    for i=1:Nconc
        GSB_down(:,i)    = VolumeData_GSB{i}(:,3);
        GSB_up(:,i)      = VolumeData_GSB{i}(:,4);
    end
    dlmwrite([scriptdir filesep subdir 'GSB_DOWN.dat'],[[0 ConcPercent'];[t2delays GSB_down.*decay]]);
    dlmwrite([scriptdir filesep subdir 'GSB_UP.dat'],[[0 ConcPercent'];[t2delays GSB_up.*decay]]);
    for i=1:Nconc
        GSB_diag12(:,i)  = VolumeData_GSB{i}(:,2);
        GSB_diag12(:,i)  = VolumeData_GSB{i}(:,1);
    end
    dlmwrite([scriptdir filesep subdir 'GSB_diag12.dat'],[[0 ConcPercent'];[t2delays GSB_diag12.*decay]]);
    dlmwrite([scriptdir filesep subdir 'GSB_diag13.dat'],[[0 ConcPercent'];[t2delays GSB_diag12.*decay]]);
    for i=1:Nconc
        ESA_diag12(:,i)  = VolumeData_ESA{i}(:,2);
        ESA_diag12(:,i)  = VolumeData_ESA{i}(:,1);
    end
    dlmwrite([scriptdir filesep subdir 'ESA_diag12.dat'],[[0 ConcPercent'];[t2delays ESA_diag12.*decay]]);
    dlmwrite([scriptdir filesep subdir 'ESA_diag13.dat'],[[0 ConcPercent'];[t2delays ESA_diag12.*decay]]);
end

if DoFit == 0
 return
end

%% Do kinetic model fit
fitPar0 = [3  20   70   90  0.3];
fitpar = zeros(length(fitPar0),Nconc);
for i=1:Nconc
    FW_data = [VolumeData_ESA{i}(:,1) VolumeData_ESA{i}(:,3)/10];
    BW_data = [VolumeData_ESA{i}(:,2) VolumeData_ESA{i}(:,4)/10];
    time = delays{i};

    LB      = [2  0    0    0   0  ];
    UB      = [5  50   500  500 1  ];

    ShowOutput = 0;
    DoFit = 1;
    tic
    fitpar(:,i) = Fit_EnT(time,FW_data,BW_data,fitPar0,LB,UB,ShowOutput,DoFit,fh);
    disp(['Fitting done for ' num2str(i) ' of ' num2str(Nconc) ' (' num2str(toc) 's)']);
end

fitpar = fitpar';