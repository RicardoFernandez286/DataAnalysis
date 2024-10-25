function [figNum,k,ydata] = plot_TransientSpectra(DataStruct,DET,data,overlay,figNum,Normalise,IncludeSSspec,varargin)

%% Hardcoded settings
FontSize = 18; % 18 for papers (single column figure), 14 for presentations (can vary)
ResPlotStyle = 'o';

switch DataStruct.Xunits
    case 'nm'
        MkrFaceAlpha = 0.15;
    case {'cm^{-1}','cm-1'}
        MkrFaceAlpha = 0.4;
    otherwise
        MkrFaceAlpha = 0.25;
end

%% Decide whether we want residuals or not
if ~isempty(varargin)
    resplot = varargin{1};
    if resplot == 0
        overlay = 1;
    end
end

%% Get the data to be plotted
% Get the TIME values to be plotted
value = DataStruct.SelTraces(:,1);

% k = index of the closest match for each selected trace
k = findClosestId2Val(DataStruct.delays,transpose(value));
k = unique(k,'stable');
L = length(k);
ydata   = zeros(size(data,1),L);
caption = strings(L,1);

for n=1:L
    % Check the delays and convert the labels and numbers to appropriate units
    if DataStruct.delays(k(n)) >= 1000 && ~strcmp(DataStruct.timescale,'s')
       newdelays = DataStruct.delays(k(n))/1000;
       switch DataStruct.timescale
           case 'fs'
               newtimescale = 'ps';
           case 'ps'
               newtimescale = 'ns';
           case 'ns'
               newtimescale = ['\mu' 's'];
       end
    elseif DataStruct.delays(k(n)) < 1 && ~strcmp(DataStruct.timescale,'s')
       newdelays = DataStruct.delays(k(n)).*1000;
       switch DataStruct.timescale
           case 'ps'
               newtimescale = 'fs';
           case 'ns'
               newtimescale = 'ps';
       end
    else
       newdelays    = DataStruct.delays(k(n));
       newtimescale = DataStruct.timescale;
    end

    % Get the actual vectors from the Z matrix (one by one)
    % together with their legend captions and corresponding ylabel
    switch Normalise
        case 0
            ydata(:,n) = data(:,k(n));
            label = '\DeltaAbs (mOD)';
        case 1
            ydata(:,n) = data(:,k(n))./max(max(abs(data(:,k(n)))));
            label = 'Normalised \DeltaAbs';
    end
    caption(n,1) = string([num2str(round(newdelays,2,'significant')) ' ' newtimescale]);
end

switch IncludeSSspec
    case 1
        % Create a new figure with consistent format
        try
            scale   = str2double(readParam('SSpecPlotSizePercent',app.SettingsPath))/100;
            filledSS= logical(str2double(readParam('SSpecPlotFillArea',app.SettingsPath)));
            faceCol = str2num(readParam('SSpecPlotFillAreaCol',app.SettingsPath));
        catch
            scale   = 60/100;
            filledSS= false;
            faceCol = 0.9*[1 1 1];
        end
        fh = figure();
        fh.Units        = 'pixels';
        fh.Position(2)  = fh.Position(2)-450*scale;
        fh.Position(3)  = 890;
        fh.Position(4)  = 425+325*scale;
        fh.Color        = [1 1 1];
        % Create the SS axis
        ax1             = axes('Parent',fh);
        % Create the transient spectra axis
        ax2             = axes('Parent',fh);
        linkaxes([ax1,ax2],'x');
        where = ax2;
    case 0
        switch overlay
            case 0
                % Create a new figure with consistent format
                fh = figure();
                fh.Name = [DataStruct.rawcorr ' DATA - TRANSIENT SPECTRA of "' DataStruct.datafilename '"'];
                fh.Position(3)  = 890;
                fh.Position(4)  = 425;
                fh.Color        = [1 1 1];
                % Define the axes
                ax2 = axes('Parent',fh);
                axes(ax2);
                ax2.Units     = 'pixels';
            case 1
                ax2=gca;
            case 2
                % Create a new figure with consistent format
                fh = figure();
                fh.Name = [DataStruct.rawcorr ' DATA - TRANSIENT SPECTRA of "' DataStruct.datafilename '"'];
                fh.Position(3)  = 890;
                fh.Position(4)  = 425;
                fh.Color        = [1 1 1];
                % Define the axes
                ax2 = axes('Parent',fh);
                axes(ax2);
                ax2.Units     = 'pixels';
        end

        where = ax2;
end

% Check that I have a fit to plot, otherwise don't do it
if overlay == 2
    if ~isfield(DataStruct,'FitData') || isempty(DataStruct.FitData{DET})
        overlay = 0;
    else
        FitData = DataStruct.FitData{DET};
        switch Normalise
            case 0
                yFit = FitData(k,:)';
            case 1
                for n=1:L
                    maxval = max(data(k(n),:));
                    minval = min(data(k(n),:));
                    if abs(maxval) >= abs(minval)
                        yFit(:,n) = FitData(k(n),:)'./maxval;
                    else
                        yFit(:,n) = FitData(k(n),:)'./minval;
                    end
                end
        end
    end
end

cm=colormap(where,othercolor('Mrainbow',L));

% Automatically detect discontinuities in the probe axis (e.g. masked pump scatter)
% and plot them accordingly
probe   = DataStruct.cmprobe{DET};
dProbe  = diff(probe) - mean(diff(probe));
jump_ID = dProbe >= 5*mean(diff(probe));
ydata(jump_ID,:) = NaN;

if overlay==2
    res(jump_ID,:)  = NaN;
    yFit(jump_ID,:) = NaN;
end

% Plot the data
hold(where,'on')
if overlay==2
    for n=1:L
       scatter(where,probe,ydata(:,n),ResPlotStyle,'MarkerEdgeColor','none','MarkerFaceColor',cm(n,:),'MarkerFaceAlpha',MkrFaceAlpha,'HandleVisibility','off');
       plot(where,probe,yFit(:,n),'LineWidth',1.5,'Marker','none','color',cm(n,:),'DisplayName',caption(n));
    end
else
    for n=1:L
       plot(where,probe,ydata(:,n),'LineWidth',1.5,'Marker','none','color',cm(n,:),'DisplayName',caption(n));
    end
end

hold(where,'off')
% Nice formatting
xlabel(where,[ DataStruct.probeunits ' (' DataStruct.Xunits ')'],'FontSize',13,'FontWeight','bold')
ylabel(where,label,'FontSize',FontSize,'FontWeight','bold')
axis(where,'tight')

if L > 15
    cb = colorbar(where);
    colormap(where,othercolor('Mrainbow',L));
    caxis(where,[min(DataStruct.delays(k)) max(DataStruct.delays(k))]);
    ylabel(cb,['Delay (' DataStruct.timescale ')'],'FontWeight','bold');
    cb.LineWidth       = 0.1;
    cb.TickLength      = 0.025;
    cb.TickDirection   = 'both';
else
    legend(where,'show');
    legend('boxoff')
    legend('Location','bestoutside')
    legend(ax2,{},'FontSize',FontSize)
end

where.Units       = 'pixels';
where.Position    = [80 75 675 320];
where.Units       = 'normalized';
where.FontSize=FontSize;

% Refline
hline = yline(where,0,'HandleVisibility','off'); hline.Color = [0.5 0.5 0.5];
ylim(where,'padded');

% Nice formatting for the combined plot
if IncludeSSspec == 1
    xl          = xlim(where);
    xl_idx      = findClosestId2Val(DataStruct.SSx,xl);
    xl_idx      = sort(xl_idx);
    SSy         = DataStruct.SSy - min(DataStruct.SSy(xl_idx(1):xl_idx(2)));

    if Normalise == 0 && max(abs(SSy(xl_idx(1):xl_idx(2)))) ~= 1
            Yname   = 'Abs. (mOD)';
    elseif Normalise == 1
            SSy     = SSy./max(abs(SSy(xl_idx(1):xl_idx(2))));
            Yname   = 'Norm. Abs.';
    elseif max(abs(SSy(xl_idx(1):xl_idx(2)))) == 1
            Yname   = 'Norm. Abs.';
    end

    minmaxIRy   = [min(SSy(xl_idx(1):xl_idx(2))) max(SSy(xl_idx(1):xl_idx(2)))];
    
    % Convert to OD if more than 500 mOD
    if max(abs(minmaxIRy)) >= 500
        SSy     = SSy./1000;
        Yname   = 'Abs. (OD)';
    end

    switch DataStruct.probeunits
        case 'Wavelength'
            SSpecName = 'UV-Vis';
        case 'Wavenumbers'
            SSpecName = 'FT-IR';
        otherwise
            SSpecName = 'SS';
    end
    
    if filledSS == 0
        plot(ax1,DataStruct.SSx,SSy,'k','LineWidth',1.5,'DisplayName',SSpecName);
    else
        area(ax1,DataStruct.SSx,SSy,'LineWidth',1.5,'EdgeColor','k','FaceColor',faceCol,'DisplayName',SSpecName);
    end

    xlim(ax1,xl);
    ylim(ax1,minmaxIRy);
    set(ax1, 'XTickLabel', []); % To remove only the labels, otherwise use XTick
%   set(ax1, 'XTick', []);
    ax1.FontSize = FontSize;
    ax2.FontSize = FontSize;
    ylabel(ax1,Yname,'FontWeight','bold','FontSize',FontSize)
    legend(ax1,'show');
    legend(ax1,'boxoff')
    legend(ax1,'Location','bestoutside')
    hline2 = yline(ax1,0,'HandleVisibility','off'); hline2.Color = [0.5 0.5 0.5];
    ylim(ax1,'padded');
    ax1.Units       = 'pixels';
    ax1.Position    = [80 400 675 320*scale];
    ax1.Units       = 'normalized';
    ax1.TickLength  = [0.015 0.025];
end

where.TickLength=[0.015 0.025];
where.Box = 'on';