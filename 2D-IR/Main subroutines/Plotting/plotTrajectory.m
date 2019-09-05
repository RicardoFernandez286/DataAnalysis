function plotTrajectory(simData,plotPie,WhichTraject,DistPlot)
%% READ from simData
version				= simData.version;
trajectInfo			= simData.trajectInfo;
trajectData_2D 		= simData.trajectData_2D;
trajectData_3D 		= simData.trajectData_3D;
Nmolecules			= simData.Nmolecules;
Nsamples			= simData.Nsamples;
Radius_LJ_Re		= simData.Radius_LJ_Re;
Radius_LJ_CN		= simData.Radius_LJ_CN;
Rbox				= simData.Rbox;

%% Plotting options
% Colors
colRe12     = brighten([0.5 0.75 1],0.5);
colRe13     = brighten([1 0.75 0.5],0.5);
colRe18     = brighten([0.75 1 0.5],0.5);
colCNBz     = 0.75.*[1 1 1];
arrowcolor  = 'k';

% Isotope shifts
iso13_shift = 48;
iso18_shift = 92;

% Dipole moment magnitude (!)
mu          = 0.1;

% Other settings
cEdgeMult   = 0.5;      % Circle edge color multiplier (must be <=1)
cAlpha      = 1;      % Circle alpha (0 = transparent, 1 = opaque)
cLineW      = 1;        % Circle line width
% DistPlot    = 'All';    % Isotopte distribution pie chart: 'Single' or 'All' trajectories

%% Create figure
fh          = figure(2);
clf(fh);
ax          = axes(fh);
ax.Visible  = 'Off';
fh.Color    = [1 1 1];
fh.Position(2) = fh.Position(2) + (fh.Position(4) - fh.Position(3));
fh.Position(4) = fh.Position(3);

% Make uniform, consistent format
ax.FontSize     = 12;
ax.LineWidth    = 1;
ax.TickLength   = [0.015 0.035];
ax.Visible      = 'On';

%% Plot the trajectory
    while WhichTraject == 0
        WhichTraject = round(Nsamples*rand);
    end
    traject = squeeze(trajectData_3D(:,:,WhichTraject));
    rLJ_Re  = ones(length(traject),1)*Radius_LJ_Re;
    if version==2
        rLJ_CN  = ones(length(traject),1)*Radius_LJ_CN;
    end

    if version == 2
        iso12   = traject(traject(:,6) == 0 & traject(:,7) ~= trajectInfo(3),:);
    else
        iso12   = traject(traject(:,6) == 0,:);
    end
    iso13   = traject(traject(:,6) == -iso13_shift,:);
    iso18   = traject(traject(:,6) == -iso18_shift,:); 

    hold(ax,'on');

    if ~isempty(iso12)
        circles(ax,iso12(:,1),iso12(:,2),rLJ_Re(1:size(iso12,1)),'facecolor',colRe12,'facealpha',cAlpha,'EdgeColor',cEdgeMult*colRe12,'LineWidth',cLineW);
    end
    if ~isempty(iso13)
        circles(ax,iso13(:,1),iso13(:,2),rLJ_Re(1:size(iso13,1)),'facecolor',colRe13,'facealpha',cAlpha,'EdgeColor',cEdgeMult*colRe13,'LineWidth',cLineW);
    end
    if ~isempty(iso18)
        circles(ax,iso18(:,1),iso18(:,2),rLJ_Re(1:size(iso18,1)),'facecolor',colRe18,'facealpha',cAlpha,'EdgeColor',cEdgeMult*colRe18,'LineWidth',cLineW);
    end

    switch version
        case 2
            CNdil   = traject(traject(:,7) == trajectInfo(3),:);
            if ~isempty(CNdil)
                circles(ax,CNdil(:,1),CNdil(:,2),rLJ_CN(1:size(CNdil,1)),'facecolor',colCNBz,'facealpha',cAlpha,'EdgeColor',cEdgeMult*colCNBz,'LineWidth',cLineW);
            end
            
            arrowpoints = traject(traject(:,7) ~= trajectInfo(3),:);
            for i=1:size(arrowpoints,1)
                    plot_arrow(ax,arrowpoints(i,1),arrowpoints(i,2),arrowpoints(i,1)+arrowpoints(i,3)*Radius_LJ_Re/mu,arrowpoints(i,2)+arrowpoints(i,4)*Radius_LJ_Re/mu,'linewidth',1,'headwidth',1,'headheight',1,'color',arrowcolor,'facecolor',arrowcolor,'edgecolor',arrowcolor);
            end
        case 1
            CNdil = [];
            for i = 1:Nmolecules
                plot_arrow(ax,traject(i,1),traject(i,2),traject(i,1)+traject(i,3)*Radius_LJ_Re/mu,traject(i,2)+traject(i,4)*Radius_LJ_Re/mu,'linewidth',1,'headwidth',1,'headheight',1,'color',arrowcolor,'facecolor',arrowcolor,'edgecolor',arrowcolor);
            end
    end

    hold(ax,'off');

    xlim(ax,[0 Rbox]);
    ylim(ax,[0 Rbox]);
    box(ax,'on');
    ax.PlotBoxAspectRatio   = [1 1 1];
    ax.FontSize             = 16;

    title(ax,['Simulated surface, trajectory # ' num2str(WhichTraject,'%i') ' of ' num2str(Nsamples,'%i')]);
    xlabel(ax,['x (' char(197) ')'],'FontWeight','bold');
    ylabel(ax,['y (' char(197) ')'],'FontWeight','bold');

%% Plot the pie chart
if strcmp(plotPie,'with pie')
    % Create figure
    fh3             = figure(3);
    clf(fh3);
    ax              = axes(fh3);
    ax.Visible      = 'Off';
    fh3.Color       = [1 1 1];
    fh3.Position(1) = fh.Position(1) + 600;
    fh3.Position(2) = fh3.Position(2) + (fh3.Position(4) - fh3.Position(3));
    fh3.Position(4) = fh3.Position(3);

    % Make uniform, consistent format
    ax.FontSize     = 12;
    ax.LineWidth    = 1;
    ax.TickLength   = [0.015 0.035];
    ax.Visible      = 'On';
    if ~isempty(CNdil)
        switch DistPlot
            case 'Single'
                MolDist     = [length(iso12) length(iso13) length(CNdil)];
            case 'All'
                MolDist     = [length(trajectData_2D(trajectData_2D(:,6) == 0,:))-length(trajectData_2D(trajectData_2D(:,7) == trajectInfo(3),:)) length(trajectData_2D(trajectData_2D(:,6) == -iso13_shift,:)) length(trajectData_2D(trajectData_2D(:,7) == trajectInfo(3),:))];
        end
    else
        switch DistPlot
            case 'Single'
                MolDist     = [length(iso12) length(iso13) length(iso18)];
            case 'All'
                MolDist     = [length(trajectData_2D(trajectData_2D(:,6) == 0,:)) length(trajectData_2D(trajectData_2D(:,6) == -iso13_shift,:)) length(trajectData_2D(trajectData_2D(:,6) == -iso18_shift,:))];
        end
    end

    pObj    = pie(ax,MolDist);
    pText   = findobj(pObj,'Type','text');
    prcntVal= get(pText,'String')';

    if ~isempty(iso18)
        colormap(ax,[colRe12;colRe13;colRe18])
        isolabel = {['Re(^{12}CO)' newline],['Re(^{13}CO)' newline],['Re(^{13}C^{18}O)' newline]};
    elseif ~isempty(CNdil)
        colormap(ax,[colRe12;colRe13;colCNBz])
        isolabel = {['Re(^{12}CO)' newline],['Re(^{13}CO)' newline],['CNBz' newline]};
    else
        colormap(ax,[colRe12;colRe13])
        isolabel = {['Re(^{12}CO)' newline],['Re(^{13}CO)' newline]};
    end

    newLabel = strcat(isolabel,prcntVal)';
    for i=1:length(pText)
        pText(i).String     = newLabel{i};
        pText(i).FontSize   = 14;
%            pText(i).FontWeight = 'bold';
    end
end