function varargout = SurfPlot_2DIR(plotOptions,dataStruct,plotaxis)

% Description: This function will make a 3D surface plot of phased 2D-IR data.
% Usage: dataStruct = ContourPlot_2DIR(plotOptions,dataStruct,plotaxis)
% Inputs:
%     plotaxis : axis handle where the plot will be made
%     
% Outputs:
%     (plot)
%     plotLimits.min/max_cut
% Ricardo Fernández-Terán / 14.04.2019 / v4.0a

debug=0;
cla(plotaxis);
%% DEBUG/MANUAL:
if debug==1
    Ncontours=50;
    ContourLineStyle='none';
    
    plot_Nwhites=0;
    plot_colourrange=100;
    plot_contourfill=1;
    plot_showcontours=1;
    
    popdelay=1;
    plot_pumpdirection='Horizontal';
    plot_axislegend=1;
    plot_colourscheme='Red/blue'; % 'Red/blue' 'Jet' 'Blue/DRed'
    plot_skiplevels=0;
    t2_delay_ps         = dataStruct.t2_delay_ps;
end

%% READ from dataStruct
% Hardcoded values
    cut_threshold       = 15; % Percentage of max intensity to cut for plotting
    LineWidth           = 0.5;
    plot_limittype      = 'Local'; % 'Global' will take min/max of the whole 2D set, while 'Local' will take only the selected region
    interpolate         = 0;
    textcolor           = 'none';
% Read data
    ProbeAxis           = dataStruct.ProbeAxis;
    PumpAxis            = dataStruct.PumpAxis;
    PROC_2D_DATA        = dataStruct.PROC_2D_DATA;
    PumpSpectrum        = dataStruct.phased_FFTZPint;
    
if debug==0
  % Read from GUI
    cut                 = plotOptions.cut;
    minWL_pump          = plotOptions.minWL_pump;
    maxWL_pump          = plotOptions.maxWL_pump;
    minWL_probe         = plotOptions.minWL_probe;
    maxWL_probe         = plotOptions.maxWL_probe;
    plot_pumpdirection  = plotOptions.plot_pumpdirection;
    plot_contourfill    = plotOptions.plot_contourfill;
    plot_Nwhites        = plotOptions.plot_Nwhites;
    Ncontours           = plotOptions.Ncontours;
    ContourLineStyle    = plotOptions.ContourLineStyle; 
    plot_axislegend     = plotOptions.plot_axislegend; % "Omega", "Pump-probe" or "w pump-probe"
    plot_colourrange    = plotOptions.plot_colourrange; % Percent of minabs/maxabs for the color range
    plot_colourscheme   = plotOptions.plot_colourscheme;
    plot_showcontours   = plotOptions.plot_showcontours;
    plot_skiplevels     = plotOptions.plot_skiplevels;
    symcolrange         = plotOptions.symcolrange;
    
    popdelay            = plotOptions.popdelay;
    
    t2delays            = dataStruct.t2delays;
    Ndelays             = length(t2delays);
    if isfield(dataStruct,'t2_startFit') ~=0
        t2_startFit     = dataStruct.t2_startFit;
    end
    
    % Read the selection
    m = popdelay;
    k = 1;

    % Determine if spectral diffusion analysis have been performed or not and whether to plot them
    if plotOptions.ShowSpecDiff && dataStruct.SpecDiff
        Plot_SpecDiff   = 1;
        % Get the spectral diffusion data from dataStruct
        CLS_Xdata   = dataStruct.CLS_Xdata;
        CLS_Ydata   = dataStruct.CLS_Ydata;
        IvCLS_Xdata = dataStruct.IvCLS_Xdata;
        IvCLS_Ydata = dataStruct.IvCLS_Ydata;
        NLS_Xdata   = dataStruct.NLS_Xdata;
        NLS_Ydata   = dataStruct.NLS_Ydata;
    else
        Plot_SpecDiff   = 0;
    end
    
    plotFitResults      = plotOptions.plotFitResults;
    plotResidualsFit    = plotOptions.plotResidualsFit;
    
    % Determine whether the data is transient 2D or not, then read the delays
    if dataStruct.Transient2D == 1
        Transient2D         = 1;
        delays              = dataStruct.t2delays(popdelay,:);
        t2_delay_ps         = delays(1);
        UV_delay_ps         = delays(2);
    else
        Transient2D         = 0;
        t2_delay_ps         = dataStruct.t2delays(popdelay,1);
    end
end

%% Plot the data
% Hide the axes for now, until everything is ready
hold(plotaxis,'off');

% Cut the datasets in half (only half the frequencies are real)
% L                   = round(length(PumpAxis{m,k})/2);
% L                   = length(PumpAxis{m,k});
% PROC_2D_DATA{m,k}   = PROC_2D_DATA{m,k}(1:L,:);
PumpAxis{m,k}       = PumpAxis{m,k}(1:L,:);
PumpSpectrum        = PumpSpectrum{m,k}(1:L,:);

% Cut the datasets along the pump direction (if selected)
if cut == 1
    cut_method = 'Axis';
else
    cut_method = 'Intensity';
end

switch cut_method
    case 'Intensity'
        accepted_pts    = PumpSpectrum >= (cut_threshold/100).*max(PumpSpectrum);
        minindex        = find(accepted_pts,1,'first');
        maxindex        = find(accepted_pts,1,'last');
    case 'Axis'
        minindex        = findClosestId2Val(PumpAxis{m,k},min(ProbeAxis))-1;
        maxindex        = findClosestId2Val(PumpAxis{m,k},max(ProbeAxis))+1;
end

switch plot_limittype
    case 'Local'
        maxabs=max(max(PROC_2D_DATA{m,k}(minindex:maxindex,:)));
        minabs=min(min(PROC_2D_DATA{m,k}(minindex:maxindex,:)));
    case 'Global'
        maxabs=max(max(PROC_2D_DATA{m,k}));
        minabs=min(min(PROC_2D_DATA{m,k}));
end

% Trim the data
proc2Ddata          = PROC_2D_DATA{m,k}(minindex:maxindex,:);
PumpAxis{m,k}       = PumpAxis{m,k}(minindex:maxindex,:);

% Get overall maxmimum of positive and negative deltaAbs
maxDabs = max([abs(minabs),abs(maxabs)]);

% % Clip the data to a limit range (i.e. set to cut value if Z â‚¬ )max_cut,min_cut( )
% if clip == 1
%     PROC_2D_DATA{m,k}(PROC_2D_DATA{m,k} > max_cut) = max_cut;
%     PROC_2D_DATA{m,k}(PROC_2D_DATA{m,k} < min_cut) = min_cut;
% end

% Cut the intensities
if symcolrange == 1
    max_cut = maxDabs*plot_colourrange/100;
    min_cut = -maxDabs*plot_colourrange/100;
else
    max_cut = maxabs*plot_colourrange/100;
    min_cut = minabs*plot_colourrange/100;
end

% Define the contours to plot
plot_contours = linspace(min_cut,max_cut,Ncontours);

% Set pump and probe orientations (i.e. define X,Y and Z)
switch plot_pumpdirection    
case 'Vertical' % Old way
        % Save XYZ
        X = ProbeAxis;
        Y = PumpAxis{m,k};
        Z = proc2Ddata;
        % Interpolate the data
        if interpolate == 1
            %X  = linspace() etcetera NOT DONE PENDING
            Z  = interp2(X,Y,Z,Xi,Yi,'bicubic');
        end
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
        Z = transpose(proc2Ddata);
        % Interpolate the data
        if interpolate == 1
            %X  = linspace() etcetera NOT DONE PENDING
            Z  = interp2(X,Y,Z,Xi,Yi,'bicubic');
        end
      % Save axes labels
        omegaX = '\omega_{1} (cm^{-1})';
        omegaY = '\omega_{3} (cm^{-1})';
        PP_x = 'Pump frequency (cm^{-1})';
        PP_y = 'Probe frequency (cm^{-1})';
        omegaPP_x = '\omega_{pump} (cm^{-1})';
        omegaPP_y = '\omega_{probe} (cm^{-1})';
end

% Make the contour plot

switch plot_contourfill
    case 1
        surfc(plotaxis,X,Y,Z,'EdgeColor','interp','FaceColor','interp','LineStyle',ContourLineStyle,'LineWidth',LineWidth);
    case 0
        surf(plotaxis,X,Y,Z,'EdgeColor','interp','FaceColor','interp','LineStyle',ContourLineStyle,'LineWidth',LineWidth);
end

% dlmwrite([dataStruct.datafilename '_traces.dat'],[[0,X'];[Y,Z]])

%% Make the plot format nice

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
colormap(plotaxis,cmap);

%% Set the colour range
% Set shading and color axis limits
shading(plotaxis,'flat');
if symcolrange == 1
    caxis(plotaxis,[-maxDabs*plot_colourrange/100,maxDabs*plot_colourrange/100]);
else
    caxis(plotaxis,[minabs*plot_colourrange/100,maxabs*plot_colourrange/100]);
end

% Show the colorbar
hcb     =   colorbar(plotaxis);
ylabel(hcb,'2D signal (a.u.)','FontWeight','bold')

% Show the zero plane on the 2D-IR spectrum
    hold(plotaxis,'on')
    surf(plotaxis,X,Y,zeros(length(Y),length(X)),'EdgeColor','none','FaceColor',0.8*[1 1 1],'FaceAlpha',0.2);
    hold(plotaxis,'off')
    
% Set appropriate limits
switch plot_pumpdirection
case 'Vertical' 
    xlim(plotaxis,[minWL_probe,maxWL_probe]);
    ylim(plotaxis,[minWL_pump,maxWL_pump]);
case 'Horizontal'
    ylim(plotaxis,[minWL_probe,maxWL_probe]);
    xlim(plotaxis,[minWL_pump,maxWL_pump]);
end

% Show the axis legends
switch plot_axislegend
    case '1'
        xlabel(plotaxis,omegaX,'FontWeight','bold','FontSize',14);
        ylabel(plotaxis,omegaY,'FontWeight','bold','FontSize',14);
    case '2'
        xlabel(plotaxis,PP_x,'FontWeight','bold','FontSize',14);
        ylabel(plotaxis,PP_y,'FontWeight','bold','FontSize',14);
    case '3'
        xlabel(plotaxis,omegaPP_x,'FontWeight','bold','FontSize',14);
        ylabel(plotaxis,omegaPP_y,'FontWeight','bold','FontSize',14);
end

% Show population delay text in upper-left corner
if abs(t2_delay_ps) < 1
    timescale   = 'fs';
    t2_delay_ps = t2_delay_ps*1000;
elseif abs(t2_delay_ps) > 1000
    timescale   = 'ns';
    t2_delay_ps = t2_delay_ps/1000;
else
    timescale   = 'ps';
end
t2_delay_ps = num2str(t2_delay_ps,'%.3g');

if Transient2D
    if abs(UV_delay_ps) < 1
    UVtimescale   = 'fs';
    UV_delay_ps = UV_delay_ps*1000;
    elseif abs(UV_delay_ps) > 1000
        UVtimescale   = 'ns';
        UV_delay_ps = UV_delay_ps/1000;
    else
        UVtimescale   = 'ps';
    end
    UV_delay_ps = num2str(UV_delay_ps,'%.3g');
    
    text(plotaxis,0.05,0.925,['t_{2} = ' t2_delay_ps ' ' timescale '; t_{UV} = ' UV_delay_ps ' ' UVtimescale],...
        'Units','normalized','FontSize',12,'FontWeight','bold','BackgroundColor',textcolor,'FontSize',14);

% If NOT transient 2D
else
    text(plotaxis,0.05,0.925,['t_{2} = ' t2_delay_ps ' ' timescale],...
        'Units','normalized','FontSize',12,'FontWeight','bold','BackgroundColor',textcolor,'FontSize',14);
end

%% Show the level lines every nth level, starting from zero
if plot_showcontours == 1
    % Define the contours to plot
    step = ((max_cut - min_cut)/Ncontours)*plot_skiplevels;
	neg_zero = (plot_Nwhites/Ncontours)*min_cut;
    pos_zero = (plot_Nwhites/Ncontours)*max_cut;

    contourlines = [min_cut:step:neg_zero pos_zero:step:max_cut];

    LineColor       = 'k';
    hold(plotaxis,'on')
    contour(plotaxis,X,Y,Z,contourlines,'LineColor',LineColor,'LineStyle',ContourLineStyle,'LineWidth',0.1);
    hold(plotaxis,'off')
end

%% Show the spectral diffusion lines
if Plot_SpecDiff
    hold(plotaxis,'on')
    % CLS
    if ~isempty(CLS_Xdata) && ~isempty(CLS_Xdata{m})
        plot(plotaxis,CLS_Xdata{m},CLS_Ydata{m}(:,popdelay),'LineWidth',2,'Color','w')
    end
    % IvCLS
    if ~isempty(IvCLS_Xdata) && ~isempty(IvCLS_Xdata{m})
        plot(plotaxis,IvCLS_Xdata{m}(:,popdelay),IvCLS_Ydata{m},'LineWidth',2,'Color','y')
    end
    % NLS
    if ~isempty(NLS_Xdata) && ~isempty(NLS_Xdata{m})
        plot(plotaxis,NLS_Xdata{m},NLS_Ydata{m},'LineWidth',2,'Color','m')
    end
    hold(plotaxis,'off')
end

% %% Show the Gaussian fit results
% if ~isempty(t2_startFit)
%     OrigDelays      = t2delays;
%     t2delays        = t2delays(t2delays>t2_startFit);
%     newDelayNumber  = popdelay-(Ndelays-length(t2delays));
%     inputstruct     = dataStruct.FitInput;
% end
% 
% if isfield(dataStruct,'FitResults') ~= 0  && ~isempty(dataStruct.FitResults) && plotFitResults && newDelayNumber > 0 
%     switch plot_pumpdirection    
%     case 'Vertical'
%         Xfit            = inputstruct.Omega{2};
%         Yfit            = inputstruct.Omega{1};
%         Zfit            = dataStruct.FitResults(:,:,newDelayNumber);
%     case 'Horizontal'
%         Xfit            = inputstruct.Omega{1};
%         Yfit            = inputstruct.Omega{2};
%         Zfit            = dataStruct.FitResults(:,:,newDelayNumber)';
%     end
%     
%     hold(plotaxis,'on')
%     contour(plotaxis,Xfit,Yfit,Zfit,plot_contours,'LineColor',0.7*[1 1 1],'LineWidth',0.1)
%     hold(plotaxis,'off')
%     
%     if plotResidualsFit == 1
%     PROC_2D_DATA    = dataStruct.PROC_2D_DATA;
%     ProbeAxis       = dataStruct.ProbeAxis;
%     PumpAxis        = dataStruct.PumpAxis{1,1};
%     % Get the indices
%     pump_range      = [min(inputstruct.Omega{1}) max(inputstruct.Omega{1})];
%     probe_range     = [min(inputstruct.Omega{2}) max(inputstruct.Omega{2})];
%     pump_idxrange   = sort(findClosestId2Val(PumpAxis,pump_range));
%     probe_idxrange  = sort(findClosestId2Val(ProbeAxis,probe_range));
%     Ndelays         = length(t2delays);
%     ZData           = zeros(pump_idxrange(2)-pump_idxrange(1)+1,probe_idxrange(2)-probe_idxrange(1)+1,Ndelays);       
%     % Cut the data
%     NOrigDelays     = length(OrigDelays);
%     for m=1:Ndelays
%         data                = PROC_2D_DATA{m+(NOrigDelays-Ndelays),1};
%         ZData(:,:,m)        = data(pump_idxrange(1):pump_idxrange(2),probe_idxrange(1):probe_idxrange(2));
%     end
%     hold(plotaxis,'on')
%     contourf(plotaxis,Xfit,Yfit,(ZData(:,:,newDelayNumber)-Zfit')',plot_contours,'LineColor','flat')
%     diagline = refline(plotaxis,1,0);
%     diagline.Color = 'k';
%     end
% end

%% Show the plot
plotaxis.Visible='On';

%% Export the plot ranges (if requested)
currentPlotLimits.max_cut = max_cut;
currentPlotLimits.min_cut = min_cut;

varargout{1} = currentPlotLimits;