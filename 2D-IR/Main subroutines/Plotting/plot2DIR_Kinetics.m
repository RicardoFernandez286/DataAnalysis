function plot2DIR_Kinetics(app,dataStruct)
% Get some information from the GUI(app structure)
plot_pumpdirection  = app.I2D_PumpAxisOrientation.Value;
interactivemode     = app.I2D_InteractivemodeSwitch.Value;
normalise           = app.I2D_NormaliseCheckBox.Value;
saveTraces          = app.I2D_SavetracesSwitch.Value;

% Get the data from dataStruct
ProbeAxis           = dataStruct.ProbeAxis;
PumpAxis            = dataStruct.PumpAxis;
Ndelays             = dataStruct.Ndelays;
t2delays            = dataStruct.t2delays;
PROC_2D_DATA        = dataStruct.PROC_2D_DATA;

% Get the desired values to plot from user
    dataStruct.SelTraces = [];
    switch interactivemode
        case 'On'
            % Select points interactively, to finish hit RETURN
            dataStruct = SelectTraces(dataStruct,0);
            if isempty(dataStruct.SelTraces)
                return
            end
            EnT = 0;
        case 'Off'
            dataStruct.SelTraces = inputdlg('Enter the coordinates of the points to plot along t2 (format: x1,y1;x2,y2...):',...
                 'Input desired coordinates', [1 60]);
            if isempty(dataStruct.SelTraces)
                return
            end
            % Plot special cases
            if strcmp(dataStruct.SelTraces,'EnT-GSB')
                dataStruct.SelTraces = [1979,1979;1979,2028;2028,2028;2028,1979];
                EnT = 1; normalise = 0;
            elseif strcmp(dataStruct.SelTraces,'EnT-ESA')
                dataStruct.SelTraces = [1979,1967;1979,2017;2028,2017;2028,1967];
                EnT = 1; normalise = 0;
            elseif strcmp(dataStruct.SelTraces,'Andrea')
                dataStruct.SelTraces = [1985,1985;1985,2062;2062,2062;2062,1985];
                EnT = 2; normalise = 0;
            else
                dataStruct.SelTraces = str2num(dataStruct.SelTraces{:}); %#ok<*ST2NM>
                EnT = 0;
            end
    end
    L = size(dataStruct.SelTraces,1);
    
% Separate the values into pump and probe, according to how the data is plotted    
switch plot_pumpdirection
    case 'Horizontal'
        pump_search     = dataStruct.SelTraces(:,1); 
        probe_search    = dataStruct.SelTraces(:,2);
    case 'Vertical'
        pump_search     = dataStruct.SelTraces(:,2); 
        probe_search    = dataStruct.SelTraces(:,1);
end

% Find the desired values in the pump and probe axis
probe_index = zeros(L,1);
pump_index  = zeros(Ndelays,L);
for i=1:L
    probe_index(i)          = findClosestId2Val(ProbeAxis,probe_search(i));
    for m=1:Ndelays
        pump_index(m,i)     = findClosestId2Val(PumpAxis{m,1},pump_search(i));
    end
end

% Initialise variables
caption = cell(L,1);
kindata = zeros(Ndelays,L);
    
for i=1:L
    % Get the values from the 2D arrays for each population delay and plot them
    % The data is by default in the format (w1,w3) (according to process2DIR.m)
    % The vector ydata will have the following format: (assuming N traces and M t2 delays)
    %         t2 delay 1: [Z1 Z2 Z3 ... Zn]
    %         t2 delay 2: [Z1 Z2 Z3 ... Zn]
    %         t2 delay 3: [Z1 Z2 Z3 ... Zn]
    %             ...
    %         t2 delay m: [Z1 Z2 Z3 ... Zn]
    
    % Get the data and normalise (if needed)
    for m=1:Ndelays
        kindata(m,i) = PROC_2D_DATA{m,1}(pump_index(m,i),probe_index(i));
    end
    
    switch normalise
        case 0
            label = '2D signal (a.u.)';
            caption{i} = ['(' num2str(round(PumpAxis{1,1}(pump_index(1,i)))) ', ' num2str(ProbeAxis(probe_index(i))) ') cm^{-1}'];
        case 1
            t0_index = findClosestId2Val(t2delays,0);
            maxval = max(kindata(t0_index:end,i));
            minval = min(kindata(t0_index:end,i));
            if abs(maxval) >= abs(minval)
                kindata(:,i) = kindata(:,i)./maxval;
                caption{i} = ['(' num2str(round(PumpAxis{1,1}(pump_index(1,i)))) ', ' num2str(ProbeAxis(probe_index(i))) ') cm^{-1}'];
            else
                kindata(:,i) = kindata(:,i)./minval;
                caption{i} = ['(' num2str(round(PumpAxis{1,1}(pump_index(1,i)))) ', ' num2str(ProbeAxis(probe_index(i))) ') cm^{-1} \times -1'];
            end
            label = 'Normalised 2D signal (a.u.)';
    end
end

% Prepare plot for EnT
if EnT ~= 0
    if abs(max(kindata(:))) <= abs(min(kindata(:)))
        kindata = -kindata;
    end
    diagFW          = kindata(:,1)./max(max(abs(kindata(:,1:2))));
    diagBW          = kindata(:,3)./max(max(abs(kindata(:,3:4))));
    xpeakFW         = 10.*kindata(:,2)./max(max(abs(kindata(:,1:2))));
    xpeakBW         = 10.*kindata(:,4)./max(max(abs(kindata(:,3:4))));
    time            = dataStruct.t2delays;
    % Create figure
    fh              = figure;
    fh.Units        = 'normalized';
    fh.Position(2)  = 0.1;
    fh.Position(4)  = fh.Position(4)*2;

    fh.Color        = [1 1 1];
    fh.Units        = 'pixels';
    ax_FW           = subplot(2,1,1);
    ax_BW           = subplot(2,1,2);
    ax_FW.FontSize  = 16;
    ax_BW.FontSize  = 16;

    box(ax_FW,'on');
    box(ax_BW,'on');

    hold(ax_FW,'on');
    hold(ax_BW,'on');
    switch EnT
        case 1
            diag1 = 'Diagonal Re(^{13}CO)';
            diag2 = 'Diagonal Re(^{12}CO)';
            xpeak1= 'Re(^{13}CO) \rightarrow Re(^{12}CO) \rm{\times10}';
            xpeak2= 'Re(^{13}CO) \leftarrow Re(^{12}CO) \rm{\times10}';
            titFW = 'Forward energy transfer';
            titBW = 'Backward energy transfer';
        case 2
            diag2 = 'Diagonal SCN^{-}';
            diag1 = 'Diagonal S^{13}C^{15}N^{-}';
            xpeak2= 'SCN^{-} \rightarrow S^{13}C^{15}N^{-} \rm{\times 10}';
            xpeak1= 'SCN^{-} \leftarrow S^{13}C^{15}N^{-} \rm{\times 10}';
            titBW = 'Downhill transfer';
            titFW = 'Uphill transfer';
    end
    % Plot the diagonal peaks
    plot(ax_FW,time,diagFW,'-or','LineWidth',1,'DisplayName',diag1);
    plot(ax_BW,time,diagBW,'-ob','LineWidth',1,'DisplayName',diag2);

    % Plot the cross peaks
    plot(ax_FW,time,xpeakFW,'-^r','LineWidth',1,'DisplayName',xpeak1);
    plot(ax_BW,time,xpeakBW,'-vb','LineWidth',1,'DisplayName',xpeak2);

    % Set axes limits
    axis(ax_FW,'tight');
    axis(ax_BW,'tight');
    
    %%% Nice formatting
    % Titles
    title(ax_FW,titFW,'FontSize',16);
    title(ax_BW,titBW,'FontSize',16);

    % Axis labels
    xlabel(ax_FW,'t_2 delay (ps)','FontWeight','bold','FontSize',16);
    xlabel(ax_BW,'t_2 delay (ps)','FontWeight','bold','FontSize',16);

    ylabel(ax_FW,'Normalised 2D signal','FontWeight','bold','FontSize',16);
    ylabel(ax_BW,'Normalised 2D signal','FontWeight','bold','FontSize',16);

    % Axis limits
    xlim(ax_FW,[0 time(end)]);
    xlim(ax_BW,[0 time(end)]);

    ylim(ax_FW,[-0.1 1.1]);
    ylim(ax_BW,[-0.1 1.1]);
    % Add zero line
    hline_FW = yline(ax_FW,0); hline_FW.Color = [0.5 0.5 0.5];
    set(get(get(hline_FW,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    hline_BW = yline(ax_BW,0); hline_BW.Color = [0.5 0.5 0.5];
    set(get(get(hline_BW,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    
    % Legends
    legend(ax_FW,'show');
    legend(ax_FW,'boxoff','FontWeight','bold')
    legend(ax_FW,{},'FontWeight','bold')

    legend(ax_BW,'show');
    legend(ax_BW,'boxoff')
    legend(ax_BW,{},'FontWeight','bold')
    % Make figure resizable
    fh.Units        = 'normalized';
else
    % Create a new figure with consistent format
    fh = figure();
    fh.Position(3)  = 800;
    fh.Position(4)  = 425;
    fh.Color        = [1 1 1];
    % Define the axes
    axes2 = axes('Parent',fh);
    axes(axes2);
    cmap=colormap(axes2,othercolor('Mrainbow',L));
    % Plot the data
    for n=1:L
       plot(axes2,dataStruct.t2delays,kindata(:,n),'LineWidth',2,'Marker','o','MarkerSize',2,'color',cmap(n,:));
       hold on
    end
    %%% Nice formatting
    axes2.FontSize = 14;
    xlabel(axes2,'t_{2} delay (ps)','FontSize',14,'FontWeight','bold');
    ylabel(axes2,label,'FontSize',14,'FontWeight','bold')
    % title([handles.datafilename;'WAITING TIME KINETICS';''],'Interpreter','none','FontSize',10)
    axis(axes2,'tight');

    % Show only positive t2 times
    xlim(axes2,[0 max(t2delays)]);

    % Create legend
    legend(axes2,caption,'FontSize',12)
    legend(axes2,'boxoff')
    legend(axes2,'Location','northeast')

    % Set linear or log X scale
    set(axes2,'xscale','lin');

    % Add zero line
    hline = yline(axes2,0); hline.Color = [0.5 0.5 0.5];
    set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

    axes2.Units     = 'pixels';
    axes2.Position  = [75 75 675 320];
    axes2.Units     = 'normalized';
end


% Decide whether to save the plotted traces or not
if strcmp(saveTraces,'On')
%     wavenumbers=transpose(dataStruct.cmprobe(k));
    if EnT == 0
        filename    = [dataStruct.rootdir filesep dataStruct.datafilename '_traces.dat'];
        data        = [[[0;PumpAxis{1,1}(pump_index(1,:))],[0;ProbeAxis(probe_index)]]';[dataStruct.t2delays(2:end),kindata(2:end,:)]];
    else
        filename    = [dataStruct.rootdir filesep dataStruct.datafilename '_EnT_traces.dat'];
        data        = [[[0;PumpAxis{1,1}(pump_index(1,:))],[0;ProbeAxis(probe_index)]]';[dataStruct.t2delays,[diagFW xpeakFW diagBW xpeakBW]]];
    end
    dlmwrite(filename,data);
end

end
