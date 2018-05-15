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
cmprobe             = handles.cmprobe;
timescale           = handles.timescale;
rawcorr             = handles.rawcorr;
plotranges          = handles.plotranges;
plot_skiplevels     = str2double(handles.plot_skiplevels.String);
plot_percentwhites  = handles.percentwhites; % Already a double, converted before
plot_showcontours   = handles.ShowContoursTick.Value;
zlimits             = str2num(handles.editMaxZ.String);
symcolrange         = handles.SymmetricColourRange_tick.Value;
LineStyle           = char(handles.ContourLineStyle.String);

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
step = (max_cut - min_cut)./Ncontours;
plot_contours = min_cut:step:max_cut;

%% Do the plot
axes(where);
X  = cmprobe;
Y  = delays;
Z  = Zdata;
contourf(where,X,Y,Z,plot_contours,'LineStyle','-','LineColor','flat');

%% Make the plot format nice:
% Colormap
n_tot       = Ncontours;
n_whites    = round(plot_percentwhites*n_tot/100);
if n_whites == 0
    n_whites = 1;
end
% If the colour range is symmetric
if symcolrange == 1
    percent_red     = 0.5-plot_percentwhites/200;
    percent_blue    = 0.5-plot_percentwhites/200;
else
    percent_red     = abs(maxabs/(abs(maxabs-minabs))) - plot_percentwhites/200;
    percent_blue    = abs(minabs/(abs(maxabs-minabs))) - plot_percentwhites/200;
end
n_blues     = round(n_tot*percent_blue);
n_reds      = round(n_tot*percent_red);
%% Create the colormaps

% Load the red map
map_red     = othercolor('YlOrRd9',n_reds);
% % Load the blue map
% map_blue    = flipud(othercolor('PuBu9',n_blues));

% Build the blue map from the red map
ylord       = othercolor('YlOrRd9',n_blues);
map_blue    = flipud([ylord(:,3) ylord(:,2) ylord(:,1)]);

% Build the white map
map_white   = ones(n_whites,3);

% Renormalize the colour scale to make the white more striking
for j=1:3
    map_blue(:,j) = map_blue(:,j)/(max(map_blue(:,j)));
    map_red(:,j) = map_red(:,j)/(max(map_red(:,j)));
end
% Join the colormap "parts" and set colormap
map = [map_blue(1:end-1,:);map_white;map_red(2:end,:)];
colormap(map);
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
xlabel(where,['Wavenumbers' ' (cm^{-1})'],'FontWeight','bold')
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