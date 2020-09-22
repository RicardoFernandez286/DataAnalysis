function plotRDF(simData,rootdir,datafilename,save)
%% READ from simData
version				= simData.version;
trajectInfo			= simData.trajectInfo;
trajectData_3D 		= simData.trajectData_3D;
Nsamples            = simData.Nsamples;
Nmolecules			= simData.Nmolecules;
Rbox				= simData.Rbox;

%% Hardcoded settings
nGR         = 200;
filt_gR     = 0;
plotArea    = 0;
plotCoordN  = 0;
plotgR_div  = 1;
% Exponent for plot
k           = 5; % can be 0 or 5

% Isotope shifts
iso13_shift = 48;
iso18_shift = 92;

%% Plot settings
% Colors
colRe12     = [0.5 0.75 1];
colRe13     = [1 0.75 0.5];
colRe18     = [0.75 1 0.5];
colCNBz     = 0.75.*[1 1 1];

%% Initialise variables
PBC                 = 1;
BoxSizes            = [Rbox, Rbox];
isoShift            = [0 iso13_shift iso18_shift];

%% Do the calculation
nTo = numel(trajectData_3D(:,1,:));
n12 = sum(squeeze(trajectData_3D(:,6,:) == -isoShift(1) & trajectData_3D(:,7,:) == trajectInfo(2)),'All');
n13 = sum(squeeze(trajectData_3D(:,6,:) == -isoShift(2) & trajectData_3D(:,7,:) == trajectInfo(2)),'All');
n18 = sum(squeeze(trajectData_3D(:,6,:) == -isoShift(3) & trajectData_3D(:,7,:) == trajectInfo(2)),'All');
nCN = sum(squeeze(trajectData_3D(:,7,:) == trajectInfo(3)),'All');

if version >= 2 && nCN > 0
    [gx,gTot,gAA,gBB,gCC,gAB,gAC,gBC] = calcPGR(trajectData_3D,nGR,BoxSizes,Nmolecules,PBC,isoShift,trajectInfo(3));
    dilName = '\otimes';
else
    [gx,gTot,gAA,gBB,gCC,gAB,gAC,gBC] = calcPGR(trajectData_3D,nGR,BoxSizes,Nmolecules,PBC,isoShift);
    dilName = '18';
end

%% Create new figure
fh                  = figure(4);
fh.Color            = 'w';
clf(fh);
axGR                = axes('parent',fh);
axGR.FontSize       = 16;
box(axGR,'on');

%% Filter and plot
if filt_gR==1
    gTot    = sgolayfilt(gTot,2,9);
    gAA     = sgolayfilt(gAA,2,9);
    gBB     = sgolayfilt(gBB,2,9);
    gCC     = sgolayfilt(gCC,2,9);
    gAB     = sgolayfilt(gAB,2,9);
    gAC     = sgolayfilt(gAC,2,9);
    gBC     = sgolayfilt(gBC,2,9);
end

hold(axGR,'on')

% plot(axGR,gx,gTot,'k','DisplayName','Total','LineWidth',2);
% plot(axGR,gx,gAA./(gx.^k),'--','Color','b','DisplayName','12-12','LineWidth',1.5);
% plot(axGR,gx,gBB./(gx.^k),'--','Color','r','DisplayName','13-13','LineWidth',1.5);
% plot(axGR,gx,gCC./(gx.^k),'--','Color','g','DisplayName',[dilName '-' dilName],'LineWidth',1.5);
plot(axGR,gx,gAB,'-','Color',colRe12.*colRe13,'DisplayName','12-13','LineWidth',2);
% plot(axGR,gx,gAC./(gx.^k),'--','Color',colRe12.*colRe18,'DisplayName',['12-' dilName],'LineWidth',1);
% plot(axGR,gx,gBC./(gx.^k),'--','Color',colRe13.*colRe18,'DisplayName',['13-' dilName],'LineWidth',1);

yline(axGR,1,'HandleVisibility','off');
% if k==0
%     yline(axGR,1,'HandleVisibility','off');
% end

if plotgR_div
    yyaxis(axGR,'right');
    plot(axGR,gx,gAB./(gx.^k),'r','DisplayName','Total','LineWidth',2);
end

% % Plot all of those with 18
% plot(axGR,gx,(gCC+gAC+gBC)./(gx.^k),'Color',colRe12.*colRe18,'DisplayName','12-18 + 13-18','LineWidth',1);
% 
% p=p+1;
% int(p) = sum(gAB./(gx.^k),1);

hold(axGR,'off')

axis(axGR,'tight')

% legend(axGR);
% legend('boxoff')

xlabel(axGR,['r* (' char(197) ')'],'FontWeight','bold');
switch k
    case 0
        ylabel(axGR,'g(r)','FontWeight','bold');
    case 5
        yyaxis(axGR,'right');
        ylabel(axGR,'g(r)/r^5','FontWeight','bold');
        yyaxis(axGR,'left');
        ylabel(axGR,'g(r)','FontWeight','bold');
end

%%
% Select which function to plot
gIJ = gAB;

% Calculate the first minimum
[~,P]   = islocalmin(gIJ,'FlatSelection','center','SamplePoints',gx);
[~,id]  = max(P);
r_min   = gx(id);
r_min   = r_min(1);

xline(axGR,r_min,'-','HandleVisibility','off','LineWidth',1,'Color',[1 0.5 0]);

endInt  = findClosestId2Val(gx,r_min);
% N_ab    = min([n12,n13])/Nsamples;
rho     = (n12+n13+n18)/Nsamples/prod(BoxSizes);

N_r = zeros(length(gx),1);
for i=2:length(gx)
 N_r(i) = 2*pi*rho*trapz(gx(1:i),gx(1:i).*gIJ(1:i));
end

if plotCoordN == 1
    fh =figure(5);
    clf(fh);
    axNr=axes('parent',figure(5));

    plot(axNr,gx,N_r,'b','LineWidth',2);
    xline(axNr,r_min,'-','HandleVisibility','off','LineWidth',1,'Color',[1 0.5 0]);
    axNr.FontSize = 16;
    xlabel(axNr,'r* (Å)','FontWeight','bold','FontSize',16);
    ylabel(axNr,'N(r)','FontWeight','bold','FontSize',16);
    % disp(max(N_r))
    axis(axNr,'tight');
    % yline(6,'--');
end
CoordNum = N_r(endInt);
disp({['R_1min = ' num2str(r_min)];['Coordination number: ' num2str(CoordNum)]})

%% Plot the area of the first coordination sphere
if plotArea == 1
    hold(axGR,'on')
    area(axGR,gx(1:endInt),gAB(1:endInt)./(gx(1:endInt).^k),'FaceColor',colRe12.*colRe13,'FaceAlpha',0.5,'HandleVisibility','off')
    hold(axGR,'off')
end

%% Save the calculated pRDF
switch save
    case 'On'
        dlmwrite([rootdir filesep datafilename '_pRDF.dat'],[gx,gTot,gAA,gBB,gCC,gAB,gAC,gBC],'delimiter','\t','precision','%.8E');
end

%% Other stuff

% disp(['Integral k=6: ' num2str((4*pi*rho)*trapz(gx(1:endInt),(gx(1:endInt).^(2-5).*(gAB(1:endInt)))))])

% sum(trajectData_2D(:,6) == 0)
% sum(trajectData_2D(:,6) == -iso13_shift)
% sum(trajectData_2D(:,6) == -iso18_shift)