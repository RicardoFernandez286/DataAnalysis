function plotRDF(simData)
%% READ from simData
version				= simData.version;
trajectInfo			= simData.trajectInfo;
trajectData_3D 		= simData.trajectData_3D;
Nmolecules			= simData.Nmolecules;
Rbox				= simData.Rbox;

%% Hardcoded settings
nGR         = 200;
filt_gR     = 0;

% Exponent for plot
k           = 0; % can be 0 or 6

% Isotope shifts
iso13_shift = 48;
iso18_shift = 92;

%% Plot settings
% Colors
colRe12     = [0.5 0.75 1];
colRe13     = [1 0.75 0.5];
colRe18     = [0.75 1 0.5];
% colCNBz     = 0.75.*[1 1 1];

%% Initialise variables
PBC                 = 1;
BoxSizes            = [Rbox, Rbox];
isoShift            = [0 iso13_shift iso18_shift];

%% Do the calculation
isCN    = squeeze(trajectData_3D(:,7,:) == trajectInfo(3));
nCN     = sum(isCN(:));

if version == 2 && nCN > 0
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

plot(axGR,gx,gTot./(gx.^k),'k','DisplayName','Total','LineWidth',2);
plot(axGR,gx,gAA./(gx.^k),'--','Color','b','DisplayName','12-12','LineWidth',1.5);
plot(axGR,gx,gBB./(gx.^k),'--','Color','r','DisplayName','13-13','LineWidth',1.5);
plot(axGR,gx,gCC./(gx.^k),'--','Color','g','DisplayName',[dilName '-' dilName],'LineWidth',1.5);
plot(axGR,gx,gAB./(gx.^k),'-','Color',colRe12.*colRe13,'DisplayName','12-13','LineWidth',2);
plot(axGR,gx,gAC./(gx.^k),'--','Color',colRe12.*colRe18,'DisplayName',['12-' dilName],'LineWidth',1);
plot(axGR,gx,gBC./(gx.^k),'--','Color',colRe13.*colRe18,'DisplayName',['13-' dilName],'LineWidth',1);

yline(axGR,0,'HandleVisibility','off');

% % Plot all of those with 18
% plot(axGR,gx,(gCC+gAC+gBC)./(gx.^k),'Color',colRe12.*colRe18,'DisplayName','12-18 + 13-18','LineWidth',1);
% 
% p=p+1;
% int(p) = sum(gAB./(gx.^k),1);

hold(axGR,'off')

axis(axGR,'tight')

legend(axGR);
legend('boxoff')

xlabel(axGR,['r* (' char(197) ')'],'FontWeight','bold');
switch k
    case 0
        ylabel(axGR,'g(r)','FontWeight','bold');
    case 6
        ylabel(axGR,'g(r)/r^6','FontWeight','bold');
end

% endInt = findClosestId2Val(gx,8.36);
% (2*pi)*trapz(gx(1:endInt),(gx(1:endInt).^2).*gAB(1:endInt))

% sum(trajectData_2D(:,6) == 0)
% sum(trajectData_2D(:,6) == -iso13_shift)
% sum(trajectData_2D(:,6) == -iso18_shift)