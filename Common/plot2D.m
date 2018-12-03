function plot2D(handles,Zdata,where,TitleOnOff)
% Description: This function will make a contour plot of transient data
% Usage: handles = plot2D(handles,Zdata,where,TitleOnOff)
% Inputs:
%     handles:      data structure 
%     Zdata:        Transient data. Rows=delays, columns=wavenumbers
%     where:        axis handle where the plot will be made
%     TitleOnOff:   switch to show (or not) the title in the plots ('On' or 'Off')
%     
% Outputs:
%     plot
% 
% Ricardo Fernández-Terán / 08.05.2018 / v2.0a

%% READ from handles
datafilename        = handles.datafilename;
delays              = handles.delays;
probe               = handles.cmprobe;
timescale           = handles.timescale;
rawcorr             = handles.rawcorr;
plotranges          = handles.plotranges;
plot_skiplevels     = str2double(handles.plot_skiplevels.String);
plot_Nwhites        = handles.percentwhites; % Already a double, converted before
plot_showcontours   = handles.ShowContoursTick.Value;
zlimits             = str2num(handles.editMaxZ.String);
symcolrange         = handles.SymmetricColourRange_tick.Value;
LineStyle           = char(handles.ContourLineStyle.String);
plot_colourscheme   = char(handles.plot_colourscheme.String{handles.plot_colourscheme.Value});

%% Parse the plot ranges
delaylim            = plotranges(1:2); %delaylim = [Min Max]
WLlim               = plotranges(3:4); %WLlim = [Min Max]
Zminmax             = plotranges(5); %Zminmax = max(|signal|) -> symmetric range
Ncontours           = plotranges(6);

% Determine the plotting range
minabs              = min(min(Zdata));
maxabs              = max(max(Zdata));
maxDabs             = max([abs(minabs),abs(maxabs)]);

if symcolrange == 1
    max_cut = maxDabs;
    min_cut = -maxDabs;
else
    max_cut = maxabs;
    min_cut = minabs;
end

% Define the contours to plot
plot_contours = linspace(min_cut,max_cut,Ncontours);

%% Do the plot
axes(where);
X  = probe;
Y  = delays;
Z  = Zdata;
contourf(where,X,Y,Z,plot_contours,'LineStyle','-','LineColor','flat');

%% Make the plot format nice:

% Create the colormaps
% Define the total number of colours and the total number of red/white/blue levels
n_tot       = Ncontours;

% Check whether the colour range is symmetric
if symcolrange == 1
    percent_red     = 0.5-plot_Nwhites/(2*Ncontours);
    percent_blue    = 0.5-plot_Nwhites/(2*Ncontours);
else
    percent_red     = abs(maxabs/(abs(maxabs-minabs)))-plot_Nwhites/(2*Ncontours);
    percent_blue    = 1 - percent_red - plot_Nwhites/(2*Ncontours);
end
n_reds      = round(n_tot*percent_red)+1;
n_blues     = round(n_tot*percent_blue)+1;
n_whites    = plot_Nwhites;

% Decide which colour scheme to use
switch plot_colourscheme
    case 'RdOr/Wh/Bl'
        % Generate the RED part of the map
        map_red     = othercolor('YlOrRd9',n_reds);
        
        % Generate the BLUE part of the map
        % This will use the same gradient as the red map, but conjugated (RGB -> BGR)
        tmp_blue    = othercolor('YlOrRd9',n_blues);
        map_blue    = flipud([tmp_blue(:,3) tmp_blue(:,2) tmp_blue(:,1)]);
%         % An alternative blue (looks like "smoked" blue)        
%         map_blue = flipud(othercolor('Blues1',n_blues));

        % Generate WHITE part of the colormap
        map_white   = ones(n_whites,3);
        
        % Renormalize the colour scale to make the white more striking
        for j=1:3
            map_blue(:,j) = map_blue(:,j)/(max(map_blue(:,j)));
            map_red(:,j) = map_red(:,j)/(max(map_red(:,j)));
        end
        % Join the colormap "parts"
        % The first/last element of the red/blue maps is white, so it is discarded
        cmap    = [map_blue(1:end-1,:);map_white;map_red(2:end,:)];
    case 'DRd/Ylw/Bl'
        cmap    = othercolor('BuDRd_12',n_tot);
    case 'Rd/Wh/Bl'
        cmap    = b2r(min_cut,max_cut,n_tot,n_whites);
    case 'DkRd/Wh/DkBl'
        cmap    = darkb2r(min_cut,max_cut,n_tot,n_whites);
    case 'Jet'
        cmap    = jet(n_tot);
end

% Set the colormap
colormap(cmap);
shading flat;

% Set color axis limits
if symcolrange == 1
    caxis([-maxDabs,maxDabs]);
else
    caxis([minabs,maxabs]);
end

% Show the colorbar
hcb=colorbar;
hcb.FontSize=12;
hcb.LineWidth = 0.1;
hcb.TickLength=0.0125;
title(hcb,{'\DeltaAbs';'(mOD)'},'FontWeight','bold','FontSize',10,'FontWeight','normal');

% Label other axes
xlabel(where,handles.probeunits,'FontWeight','bold')
ylabel(where,['Delays' ' (' timescale ')'],'FontWeight','bold')
hline = refline(0,0); hline.Color = [0.5 0.5 0.5];
if strcmp(TitleOnOff,'On')
    title({datafilename;[rawcorr,' DATA';'']},'Interpreter','none')  
end

% Adjust the limits (if given) -> MISSING "IF GIVEN" PART
xlim(WLlim)
ylim(delaylim)

% Log/linear time axis
set(where,'yscale',handles.linlog)

%% Show the level lines every nth level, starting from zero
if plot_showcontours == 1
    % Define the contours to plot
    step = ((max_cut - min_cut)./Ncontours).*plot_skiplevels;
	neg_zero = plot_percentwhites*min_cut/100;
    pos_zero = plot_percentwhites*max_cut/100;
    
    plot_contours = [min_cut:step:neg_zero pos_zero:step:max_cut];

    LineColor       = 'k';
    hold on
    contour(where,X,Y,Z,plot_contours,'LineColor',LineColor,'LineStyle',LineStyle,'LineWidth',0.1);
    hold off
end