function SurfPlot_2DIR(handles,plotaxis)

% Description: This function will plot already phased 2D-IR data.
% Usage: handles = plot_2DIR(handles,plotaxis)
% Inputs:
%     plotaxis : axis handle where the plot will be made
%     
% Outputs:
%     plot
% 
% Ricardo Fernández-Terán / 24.03.2018 / v2.1c

debug=0;
% %% DEBUG/MANUAL:
% if debug==1
%     Ncontours=50;
%     LineStyle='-';
%     LineColor='flat';
%     
%     plot_percentwhites=0;
%     plot_colorrange=100;
%     plot_contourfill=1;
%     plot_showcontours=1;
%     
%     popdelay=1;
%     plot_pumpdirection='Horizontal';
%     plot_axislegend=1;
%     plot_colourscheme='Red/blue'; % 'Red/blue' 'Jet' 'Blue/DRed'
%     plot_skiplevels=0;
%     t2_delay_ps         = handles.t2_delay_ps;
% end

%% READ from handles
% Hardcoded values
    cut                 = handles.CutPlot_tick.Value;
    LineWidth           = 0.5;
    plot_limittype      = 'Local'; % 'Global' will take min/max of the whole 2D set, while 'Local' will take only the selected region

% Read data
    ProbeAxis           = handles.ProbeAxis;
    PumpAxis            = handles.PumpAxis;
    PROC_2D_DATA        = handles.PROC_2D_DATA;
    
if debug==0
  % Read from GUI
    minWL               = str2double(handles.editWLmin.String);
    maxWL               = str2double(handles.editWLmax.String);
    plot_pumpdirection  = char(handles.plot_pumpdirection.String{handles.plot_pumpdirection.Value});
    plot_contourfill    = handles.plot_contourfill.Value;
    plot_percentwhites  = str2double(handles.plot_percentwhites.String);
    Ncontours           = str2double(handles.plot_Ncontours.String);
    LineStyle           = char(handles.plot_contourLineStyle.String); 
    plot_axislegend     = handles.plot_axislegend.Value; % "Omega", "Pump-probe" or "w pump-probe"
    plot_colorrange     = str2double(handles.plot_colorrange.String); % Percent of minabs/maxabs for the color range
    plot_colourscheme   = char(handles.plot_colourscheme.String{handles.plot_colourscheme.Value});
    plot_showcontours   = handles.ShowContoursTick.Value;
    plot_skiplevels     = str2double(handles.plot_skiplevels.String);
    popdelay            = handles.Population_delay.Value;
    symcolrange         = handles.SymColRange_tick.Value;
    
    t2_delay_ps         = num2str(str2double(char(handles.Population_delay.String(popdelay))),'%.3g');
    
    LineColor           = [1 1 1]*0;
    axes(plotaxis);
end

%% Plot the data
% Hide the axes for now, until everything is ready
plotaxis.Visible='Off';

% Read the selection
m = popdelay;
k = 1;

% Cut the datasets along the pump direction
minindex    = findClosestId2Val(PumpAxis{m,k},min(ProbeAxis))-1;
maxindex    = findClosestId2Val(PumpAxis{m,k},max(ProbeAxis))+1;
if cut == 1
    PROC_2D_DATA{m,k}   = PROC_2D_DATA{m,k}(minindex:maxindex,:);
    PumpAxis{m,k}       = PumpAxis{m,k}(minindex:maxindex,:);
end

switch plot_limittype
    case 'Local'
        if cut == 0
            % Get the LOCAL intensity min/max
            maxabs=max(max(PROC_2D_DATA{m,k}(minindex:maxindex,:)));
            minabs=min(min(PROC_2D_DATA{m,k}(minindex:maxindex,:)));
        elseif cut == 1
            % Get the GLOBAL intensity min/max (the dataset is already cut)
            maxabs=max(max(PROC_2D_DATA{m,k}));
            minabs=min(min(PROC_2D_DATA{m,k}));
        end
    case 'Global'
            maxabs=max(max(PROC_2D_DATA{m,k}));
            minabs=min(min(PROC_2D_DATA{m,k}));
end

% Get overall maxmimum of positive and negative deltaAbs
maxDabs = max([abs(minabs),abs(maxabs)]);

% Cut the intensities
if symcolrange == 1
    max_cut = maxDabs*plot_colorrange/100;
    min_cut = -maxDabs*plot_colorrange/100;
else
    max_cut = maxabs*plot_colorrange/100;
    min_cut = minabs*plot_colorrange/100;
end

% if clip == 1
%     PROC_2D_DATA{m,k}(PROC_2D_DATA{m,k} > max_cut) = max_cut;
%     PROC_2D_DATA{m,k}(PROC_2D_DATA{m,k} < min_cut) = min_cut;
% end

% Define the contours to plot
step = (max_cut - min_cut)./Ncontours;
plot_contours = min_cut:step:max_cut;

% Set pump and probe orientations (i.e. define X,Y and Z)
switch plot_pumpdirection    
case 'Vertical' % Old way
        % Save XYZ
        X = ProbeAxis;
        Y = PumpAxis{m,k};
        Z{m,k}= PROC_2D_DATA{m,k};
        % Save axes labels
        omegaX = '\omega_{3} (cm^{-1})';
        omegaY = '\omega_{1} (cm^{-1})';
        PP_x = 'Probe frequency (cm^{-1})';
        PP_y = 'Pump frequency (cm^{-1})';
        omegaPP_x = '\omega_{probe} (cm^{-1})';
        omegaPP_y = '\omega_{pump} (cm^{-1})';
case 'Horizontal' % New way
      % Save XYZ
        X = PumpAxis{m,k};
        Y = ProbeAxis;
        Z{m,k} = transpose(PROC_2D_DATA{m,k});
      % Save axes labels
        omegaX = '\omega_{1} (cm^{-1})';
        omegaY = '\omega_{3} (cm^{-1})';
        PP_x = 'Pump frequency (cm^{-1})';
        PP_y = 'Probe frequency (cm^{-1})';
        omegaPP_x = '\omega_{pump} (cm^{-1})';
        omegaPP_y = '\omega_{probe} (cm^{-1})';
end

% Make the surface plot
EdgeColor   = 'interp';
       
switch plot_contourfill
    case 1
        surfc(plotaxis,X,Y,Z{m,k},'EdgeColor',EdgeColor,'FaceColor','interp','LineStyle',LineStyle,'LineWidth',LineWidth,'FaceLighting','gouraud');
    case 0
        surf(plotaxis,X,Y,Z{m,k},'EdgeColor',EdgeColor,'FaceColor','interp','LineStyle',LineStyle,'LineWidth',LineWidth,'FaceLighting','gouraud');
end

%% Make the plot format nice

% Define how many colours
n_tot = Ncontours;

% Decide which colour scheme to use
switch plot_colourscheme
    case 'Red/blue'
        n_whites = round(plot_percentwhites*n_tot/100);
        if n_whites == 0
            n_whites = 1;
        end
        % If the colour range is symmetric
        if symcolrange == 1
            percent_red = 0.5-plot_percentwhites/200;
            percent_blue = 0.5-plot_percentwhites/200;
        else
            percent_red = abs(maxabs/(abs(maxabs-minabs)))-plot_percentwhites/200;
            percent_blue = 1 - percent_red - plot_percentwhites/200;
        end
        n_blues     = round(n_tot*percent_blue);
        n_reds      = round(n_tot*percent_red);
        % Load the colormaps
        map_red     = othercolor('YlOrRd9',n_reds);
        %map_blue = flipud(othercolor('PuBu9',n_blues));
        ylord       = othercolor('YlOrRd9',n_blues);
        map_blue    = flipud([ylord(:,3) ylord(:,2) ylord(:,1)]);
        map_white   = ones(n_whites,3);
        % Renormalize the colour scale to make the white more striking
        for j=1:3
            map_blue(:,j) = map_blue(:,j)/(max(map_blue(:,j)));
            map_red(:,j) = map_red(:,j)/(max(map_red(:,j)));
        end
        % Join the colormap "parts"
        map = [map_blue(1:end-1,:);map_white;map_red(2:end,:)];
        colormap(map);
        % Set shading and color axis limits
        shading flat;
    case 'Jet'
        colormap(jet(n_tot));
    case 'Blue/DRed'
        colormap(othercolor('BuDRd_12',n_tot));
end

% Set the colour range
if symcolrange == 1
    caxis([-maxDabs*plot_colorrange/100,maxDabs*plot_colorrange/100]);
else
    caxis([minabs*plot_colorrange/100,maxabs*plot_colorrange/100]);
end

% Show the colorbar
hcb=colorbar;
ylabel(hcb,'2D signal (a.u.)','FontWeight','bold')

% Show the diagonal line on the 2D-IR spectrum
    diagline = refline(1,0);
    diagline.Color = [0 0 0];
    diagline.LineWidth = 0.5;

% Set appropriate limits
switch plot_pumpdirection
    case 'Vertical' 
        xlim(plotaxis,[min(ProbeAxis),max(ProbeAxis)]);
        ylim(plotaxis,[minWL,maxWL]);
    case 'Horizontal'
        ylim(plotaxis,[min(ProbeAxis),max(ProbeAxis)]);
        xlim(plotaxis,[minWL,maxWL]);
end

% Show the axis legends
switch plot_axislegend
    case 1
        xlabel(omegaX,'FontWeight','bold','FontSize',14);
        ylabel(omegaY,'FontWeight','bold','FontSize',14);
    case 2
        xlabel(PP_x,'FontWeight','bold','FontSize',14);
        ylabel(PP_y,'FontWeight','bold','FontSize',14);
    case 3
        xlabel(omegaPP_x,'FontWeight','bold','FontSize',14);
        ylabel(omegaPP_y,'FontWeight','bold','FontSize',14);
end

% Show population delay text in upper-left corner
text(0.05,0.90,['t_{2} = ' t2_delay_ps ' ps'],'Units','normalized','FontSize',12); % 'FontWeight','bold'

% % Make the figure square
% pbaspect([1 1 1])
% daspect([1 1 1])

%% Show the level lines every nth level, starting from zero
if plot_showcontours == 1
    % Define the contours to plot
    step = ((max_cut - min_cut)./Ncontours).*plot_skiplevels;
	neg_zero = plot_percentwhites*min_cut/100;
    pos_zero = plot_percentwhites*max_cut/100;
    
    plot_contours = [min_cut:step:neg_zero pos_zero:step:max_cut];

    LineColor       = 'k';
    hold on
    contour(plotaxis,X,Y,Z{m,k},plot_contours,'LineColor',LineColor,'LineStyle',LineStyle,'LineWidth',0.1);
    hold off
end

% %% Show the level lines every nth level starting from zero (in + and - direction) - OLD WAY
% if plot_showcontours == 1
%     % First, find level zero
%     lvl0 = findClosestId2Val(hContour.LevelList,0);
%     drawnow
%     % Then, hide everything
%     set(hContour.EdgePrims,'Visible','Off');
%     % Next, show only the level lines that we want to plot
%     switch plot_colourscheme
%         case 'Red/blue'
%             set(hContour.EdgePrims(1:plot_skiplevels:end),'Visible','On');
%         otherwise 
%             set(hContour.EdgePrims(1:plot_skiplevels:lvl0),'Visible','On');
%             set(hContour.EdgePrims(lvl0:plot_skiplevels:end),'Visible','On');
%           % Hide explicitly the closest level to level zero
%             set(hContour.EdgePrims(lvl0),'Visible','Off');
%             % set(hContour.FacePrims(findClosestId2Val(hContour.LevelList,0)),'Visible','Off');
%     end
% end

%% Clip the limits
% % Get the new plotting limits
%     [minclip,maxclip] = caxis;
% % First make sure that everything is shown
% set(hContour.FacePrims,'Visible','On');
% set(hContour.FacePrims,'Visible','On');
% % Then, hide the faces and lines that we don't want to see
% if plot_contourfill ~= 0 && plot_colorrange < 100
%     set(hContour.FacePrims(1),'Visible','Off');
%     set(hContour.FacePrims(end),'Visible','Off');
% end
% drawnow
% set(hContour.EdgePrims(1:end),'Visible','Off');
% set(hContour.EdgePrims(1:end),'Visible','Off');

%% Show the plot
plotaxis.Visible='On';