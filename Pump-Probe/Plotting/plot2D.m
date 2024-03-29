function plot2D(app,Zdata,where,TitleOnOff,k,varargin)
% Description: This function will make a contour plot of transient data
% Usage: handles = plot2D(dataStruct,Zdata,where,TitleOnOff)
% Inputs:
%     handles:      data structure 
%     Zdata:        Transient data. Rows=delays, columns=wavenumbers
%     where:        axis handle where the plot will be made
%     TitleOnOff:   switch to show (or not) the title in the plots ('On' or 'Off')
%     
% Outputs:
%     plot
% 
% Ricardo Fernandez-Teran / 31.08.2022 / v3.1b

%% READ from Data Structure
dataStruct          = app.PP_Data;
datafilename        = dataStruct.datafilename;
delays              = dataStruct.delays;
probe               = dataStruct.cmprobe{k};
rawcorr             = dataStruct.rawcorr;
plotranges          = dataStruct.plotranges{k};

% Read plot options from the GUI app
plot_Ncontours      = app.PP_NContours.Value;
plot_skiplevels     = app.PP_NskipContours.Value;
plot_Nwhites        = app.PP_NWhites.Value;
plot_percentwhites  = plot_Nwhites/plot_Ncontours;
plot_showcontours   = app.PP_ShowContours_tick.Value;
plot_filledcontours = app.PP_FilledContours_tick.Value;
plot_scaleRange     = app.PP_PlotScaleSlider.Value/100;
symcolrange         = app.PP_SymmetricColourScale.Value;
LineStyle           = app.PP_ContourFormat.Value;
plot_colourscheme   = app.PP_ColourScheme.Value;

linlog_time         = dataStruct.linlog;

if ~isempty(varargin)
   quickPlot = varargin{1};
else
   quickPlot = false;
end

%% Parse the plot ranges
delaylim            = plotranges(1:2); %delaylim = [Min Max]
WLlim               = plotranges(3:4); %WLlim = [Min Max]

% Determine the plotting range
minabs              = min(min(Zdata{k}));
maxabs              = max(max(Zdata{k}));
maxDabs             = max([abs(minabs),abs(maxabs)]);

switch TitleOnOff
    case 'Noise'
        symcolrange = 0;
end

if symcolrange == 1
    max_cut = maxDabs*plot_scaleRange;
    min_cut = -maxDabs*plot_scaleRange;
else
    max_cut = maxabs*plot_scaleRange;
    min_cut = minabs*plot_scaleRange;
end

switch TitleOnOff
    case 'Noise'
        min_cut = 0;
end

% Define the contours to plot
plot_contours = linspace(min_cut,max_cut,plot_Ncontours);

%% Do the plot
X  = probe;
Y  = delays;
Z  = Zdata{k};

% Automatically detect discontinuities in the probe axis (e.g. masked pump scatter)
% and plot them accordingly
dProbe  = diff(probe) - mean(diff(probe));
jump_ID = dProbe >= 5*mean(diff(probe));
Z(:,jump_ID) = NaN;

if sum(abs(Zdata{k}(:)),'omitnan') == 0
    warning('Something is not right! The pump-probe data is identically zero!')
    return
end

if quickPlot || numel(Zdata{k}) >= 2048*100
    pcolor(where,X,Y,Z);
else
    if plot_filledcontours == 1
        contourf(where,X,Y,Z,plot_contours,'LineStyle','-','LineColor','flat');
    else
        contour(where,X,Y,Z,plot_contours,'LineStyle','-','LineColor','flat');
    end
end    

%% Make the plot format nice:

% Create the colormaps
% Define the total number of colours and the total number of red/white/blue levels
n_tot       = plot_Ncontours;

% Check whether the colour range is symmetric
if symcolrange == 1
    percent_red     = 0.5-plot_Nwhites/(2*plot_Ncontours);
    percent_blue    = 0.5-plot_Nwhites/(2*plot_Ncontours);
else
    percent_red     = abs(maxabs/(abs(maxabs-minabs)))-plot_Nwhites/(2*plot_Ncontours);
    percent_blue    = 1 - percent_red - plot_Nwhites/(2*plot_Ncontours);
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
    case 'Spectral'
        cmap    = flipud(othercolor('Spectral11',n_tot));
    case 'Turbo'
        cmap    = turbo(n_tot);
end

% Set the colormap
colormap(where,cmap);
shading(where,'flat');

% Set color axis limits
caxis(where,[min_cut,max_cut]);

% Show the colorbar
hcb=colorbar(where);
hcb.FontSize=12;
hcb.LineWidth = 0.1;
hcb.TickLength=0.025;
hcb.TickDirection = 'out';
title(hcb,{'\DeltaAbs';'(mOD)'},'FontWeight','bold','FontSize',10,'FontWeight','bold');

% Label other axes
% Probe_label = ['$$ \bf{' dataStruct.probeunits '~( ' dataStruct.Xunits ')}$$'];
% Time_label  = ['$$ \bf{Delays~(' dataStruct.timescale ')} $$'];
% xlabel(where,Probe_label,'Interpreter','latex');
% ylabel(where,Time_label,'Interpreter','latex');

Probe_label = [ dataStruct.probeunits ' (' dataStruct.Xunits ')'];
Time_label  = ['Delays (' dataStruct.timescale ')'];
xlabel(where,Probe_label,'FontWeight','bold','Interpreter','tex');
ylabel(where,Time_label,'FontWeight','bold','Interpreter','tex');

hline = yline(where,0); hline.Color = [0.5 0.5 0.5];
switch TitleOnOff
    case 'On'
        title(where,{datafilename;[rawcorr,' DATA';'']},'Interpreter','none')  
    case 'Noise'
        title(where,{datafilename;'NOISE PROFILE'},'Interpreter','none')
        title(hcb,{'\delta\DeltaAbs';'(mOD)'},'FontWeight','bold','FontSize',12,'FontWeight','bold');
end
% Adjust the limits (if given)
xlim(where,WLlim)
ylim(where,delaylim)

% Log/linear time axis
where.YScale = linlog_time;

%% Show the level lines every nth level, starting from zero
if plot_showcontours == 1
    % Define the contours to plot
    step = round(((max_cut - min_cut)./plot_Ncontours).*plot_skiplevels);
	neg_zero = plot_percentwhites*min_cut/100;
    pos_zero = plot_percentwhites*max_cut/100;
    
    plot_contours = [min_cut:step:neg_zero pos_zero:step:max_cut];

    LineColor       = 'k';
    hold(where,'on')
    contour(where,X,Y,Z,plot_contours,'LineColor',LineColor,'LineStyle',LineStyle,'LineWidth',0.1);
    hold(where,'off')
end

%% Show axes
where.Visible       = 'on';
where.Box           = 'on';
where.Layer         = 'top';
where.TickLength    = [0.03 0.025];
