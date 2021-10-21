function plot2DIR_Slices(app,dataStruct,k)
% Ask the user what to plot
slice_options = {'Diagonal','Along fixed pump WL','Along fixed probe WL','Integrate along pump axis','Integrate along probe axis','Along several pump WL','Along several probe WL'};
[slice_typeindx,doplot] = listdlg('ListString',slice_options,'OKstring','Plot','SelectionMode','single','ListSize',[150,120],'PromptString','Select slice type to plot:');

if doplot == 0
    return
end

% Get some information from the GUI
PopDelay            = app.I2D_t2delays.Value;
plot_pumpdirection  = app.I2D_PumpAxisOrientation.Value;
interactivemode     = app.I2D_InteractivemodeSwitch.Value;
normalise           = app.I2D_NormaliseCheckBox.Value;
cutplot             = app.I2D_CutPlot_tick.Value;
% Get the data
ProbeAxis           = dataStruct.ProbeAxis{k};
PumpAxis            = dataStruct.PumpAxis;
Ndelays             = dataStruct.Ndelays;
t2delays            = dataStruct.t2delays;
PROC_2D_DATA        = dataStruct.PROC_2D_DATA;
datafilename        = string(dataStruct.datafilename);

switch slice_options{slice_typeindx}
    case 'Diagonal'
        % Initialise variables
        pump_indexes    = cell(Ndelays,1);
        data            = zeros(length(ProbeAxis),Ndelays);
        
        % Get the data to be plotted
        for m=1:Ndelays
            % Get the indices of the probe wavelengths in the pump axis
            pump_indexes{m}     = findClosestId2Val(PumpAxis{m,1},transpose(ProbeAxis));
            probe_indexes       = 1:length(ProbeAxis);
            % Get the data
            data(:,m)           = diag(PROC_2D_DATA{m,1}(pump_indexes{m},probe_indexes));
        end
        
        % Create a new figure
        fh = figure();

        % Define the axes
        axes2 = axes(fh);
        
        % Plot the data
        cmap=colormap(axes2,othercolor('Mrainbow',Ndelays));
        for m=1:Ndelays
           plot(axes2,ProbeAxis,data(:,m),'-','LineWidth',2,'MarkerSize',2,'color',cmap(m,:),'DisplayName',[num2str(t2delays(m),'%.3g') ' ps']);
           hold(axes2,'on');
        end
        
        % Save plot details
        plot_title  = 'DIAGONAL CUTS';
        x_axis      = 'Probe';
        
        %%% Nice formatting
        % Size and colour
        fh.Color = [1 1 1];
        currsize = fh.Position;
        fh.Position = [currsize(1:2) 900 425];
        axes2.Units = 'Pixels';
        currsize = axes2.OuterPosition;
        axes2.OuterPosition = [currsize(1:2) 950 420];
        axes2.Units = 'Normalized';

        % Title, axis labels and legend
        axes2.FontSize = 14;
        xlabel(axes2,[x_axis ' wavelength (cm^{-1})'],'FontSize',14,'FontWeight','bold');
        ylabel(axes2,'2D signal (a.u.)','FontSize',14,'FontWeight','bold')
        title(axes2,[datafilename;plot_title;''],'Interpreter','none','FontSize',10)
        axis tight

        % Create legend
        legend(axes2,'show');
        legend(axes2,'boxoff')
        legend(axes2,'Location','eastoutside')

        % Add zero line
        hline = yline(axes2,0,'HandleVisibility','off');
        hline.Color = [0.5 0.5 0.5];
        
    case 'Along fixed pump WL'
        switch interactivemode
            case 'On'
                % Select points interactively, to finish hit RETURN
                dataStruct = SelectTraces(dataStruct,0); 
                % Separate the values into pump and probe, according to how the data is plotted    
                switch plot_pumpdirection
                    case 'Horizontal'
                        pump_search    = dataStruct.SelTraces(:,1); 
                    case 'Vertical'
                        pump_search    = dataStruct.SelTraces(:,2); 
                end
            case 'Off'
                dataStruct.SelTraces = inputdlg('Enter the pump wavenumbers to plot:',...
                     'Input pump wavenumbers to plot:', [1 60]);
                pump_search = str2num(dataStruct.SelTraces{:})';
        end
        L = length(pump_search);

        % Initialise variables
        pump_indexes    = cell(Ndelays,1);
        data            = zeros(length(ProbeAxis),Ndelays,L);
        
        % Get the data to be plotted
        for m=1:Ndelays
            % Get the indices of the probe wavelengths
            pump_indexes{m}     = findClosestId2Val(PumpAxis{m,1},pump_search);
            probe_indexes       = 1:length(ProbeAxis);
            % Get the data
            for p=1:L
                data(:,m,p)     = transpose(PROC_2D_DATA{m,1}(pump_indexes{m}(p),probe_indexes));
            end
        end
        
        % Plot each pump wavelength in one new figure
        for p=1:L
            % Create a new figure
            fh = figure();

            % Define the axes
            axes2 = axes('parent',fh);

            % Plot the data
            cmap=colormap(axes2,othercolor('Mrainbow',Ndelays));
            for m=1:Ndelays
               plot(axes2,ProbeAxis,data(:,m,p),'-','LineWidth',2,'MarkerSize',2,'color',cmap(m,:),'DisplayName',[num2str(t2delays(m),'%.3g') ' ps']);
               hold(axes2,'on');
            end

            % Save plot details
            plot_title  = ['TRANSIENTS AT ' num2str(PumpAxis{1,1}(pump_indexes{1}(p)),'%.4g') ' cm^{-1} (pump)'];
            x_axis      = 'Probe';
            
            %%% Nice formatting
            % Size and colour
            fh.Color = [1 1 1];
            currsize = fh.Position;
            fh.Position = [currsize(1:2) 900 425];
            axes2.Units = 'Pixels';
            currsize = axes2.OuterPosition;
            axes2.OuterPosition = [currsize(1:2) 950 420];
            axes2.Units = 'Normalized';

            % Title, axis labels and legend
            axes2.FontSize = 14;
            xlabel(axes2,[x_axis ' wavelength (cm^{-1})'],'FontSize',14,'FontWeight','bold');
            ylabel(axes2,'2D signal (a.u.)','FontSize',14,'FontWeight','bold')
            title(axes2,plot_title,'FontSize',10)
            axis(axes2,'tight');

            % Create legend
            legend(axes2,'show');
            legend('boxoff')
            legend('Location','eastoutside')

            % Add zero line
            hline = yline(axes2,0,'HandleVisibility','off');
            hline.Color = [0.5 0.5 0.5];
        end
    
    case 'Along fixed probe WL'
        % Get the desired values to plot from user
        dataStruct.SelTraces = [];
        switch interactivemode
            case 'On'
                % Select points interactively, to finish hit RETURN
                dataStruct = SelectTraces(dataStruct,0); 
                % Separate the values into pump and probe, according to how the data is plotted    
                switch plot_pumpdirection
                    case 'Horizontal'
                        probe_search    = dataStruct.SelTraces(:,2); 
                    case 'Vertical'
                        probe_search    = dataStruct.SelTraces(:,1); 
                end
            case 'Off'
                dataStruct.SelTraces = inputdlg('Enter the probe wavenumbers to plot:',...
                     'Input probe wavenumbers to plot:', [1 60]);
                probe_search = str2num(dataStruct.SelTraces{:}); %#ok<*ST2NM>
        end
        L = length(probe_search);
        
        % Get the segment of the pump axis
        pump_indexes = cell(Ndelays,1);
        for m=1:Ndelays
            if cutplot
                min_pump        = findClosestId2Val(PumpAxis{m,1},min(ProbeAxis));
                max_pump        = findClosestId2Val(PumpAxis{m,1},max(ProbeAxis));
                pump_indexes{m} = min_pump:1:max_pump;
            else
                pump_indexes{m} = 1:1:length(PumpAxis{m,1});
            end
        end
        
        % Initialise variables
        data                    = zeros(length(PumpAxis{1,1}(pump_indexes{m})),Ndelays,L);
              
        % Get the probe indexes
        probe_indexes           = findClosestId2Val(ProbeAxis,probe_search);
        
        % Get the data to be plotted
        for m=1:Ndelays
            for p=1:L
                data(:,m,p)     = PROC_2D_DATA{m,1}(pump_indexes{m},probe_indexes(p));
            end
        end
        
        % Plot each probe wavelength in one new figure
        for p=1:L
            % Create a new figure
            fh = figure();

            % Define the axes
            axes2 = axes('parent',fh);

            % Plot the data
            cmap=colormap(axes2,othercolor('Mrainbow',Ndelays));
            for m=1:Ndelays
               plot(axes2,PumpAxis{m,1}(pump_indexes{m}),data(:,m,p),'-','LineWidth',2,'MarkerSize',2,'color',cmap(m,:),'DisplayName',[num2str(t2delays(m),'%.3g') ' ps']);
               hold(axes2,'on')
            end

            % Save plot details
            plot_title  = ['TRANSIENTS AT ' num2str(ProbeAxis(probe_indexes(p)),'%.4g') ' cm^{-1} (probe)'];
            x_axis      = 'Pump';
            
            %%% Nice formatting
            % Size and colour
            fh.Color = [1 1 1];
            currsize = fh.Position;
            fh.Position = [currsize(1:2) 900 425];
            axes2.Units = 'Pixels';
            currsize = axes2.OuterPosition;
            axes2.OuterPosition = [currsize(1:2) 950 420];
            axes2.Units = 'Normalized';

            % Title, axis labels and legend
            axes2.FontSize = 14;
            xlabel(axes2,[x_axis ' wavelength (cm^{-1})'],'FontSize',12,'FontWeight','bold');
            ylabel(axes2,'2D signal (a.u.)','FontSize',12,'FontWeight','bold')
            title(axes2,plot_title,'FontSize',10)
            axis(axes2,'tight');

            % Create legend
            legend(axes2,'show');
            legend(axes2,'boxoff')
            legend(axes2,'Location','eastoutside')

            % Add zero line
            hline = yline(axes2,0,'HandleVisibility','off');
            hline.Color = [0.5 0.5 0.5];
            set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        end

    case 'Integrate along pump axis'
        switch interactivemode
            case 'Off'
                % Ask user for pump range (default is same as probe range)
                defaults  = {num2str([ProbeAxis(1) ProbeAxis(end)])};
                SelTraces = inputdlg('Enter the pump range to integrate:',...
                             'Input pump wavenumber range to integrate:', [1 60],defaults);
                SelTraces = str2num(SelTraces{:});
            case 'On'
                SelTraces = [ProbeAxis(1) ProbeAxis(end)];
        end
        
        % Get the indices of the pump wavelengths to integrate
        pump_indexes    = findClosestId2Val(PumpAxis{1,1},SelTraces);
        % Initialize variables
        data            = zeros(length(ProbeAxis),Ndelays);
        % Get the data to be plotted
        for m=1:Ndelays
            data(:,m)   = transpose(sum(PROC_2D_DATA{m,1}(min(pump_indexes):max(pump_indexes),:),1));
        end
        
        % Create a new figure
        fh = figure();

        % Define the axes
        axes2 = axes(fh);

        % Plot the data
        cmap=colormap(axes2,othercolor('Mrainbow',Ndelays));
        for m=1:Ndelays
           plot(axes2,ProbeAxis,data(:,m),'-','LineWidth',2,'MarkerSize',2,'color',cmap(m,:),'DisplayName',[num2str(t2delays(m),'%.3g') ' ps']);
           hold(axes2,'on');
        end

        % Save plot details
        plot_title  = ['INTEGRATED 2D SIGNAL FOR PUMP WL' num2str(PumpAxis{1,1}(min(pump_indexes))) ' to ' num2str(PumpAxis{1,1}(max(pump_indexes))) ' cm^{-1}'];
        x_axis      = 'Probe';

        %%% Nice formatting
        % Size and colour
        fh.Color = [1 1 1];
        currsize = fh.Position;
        fh.Position = [currsize(1:2) 900 425];
        axes2.Units = 'Pixels';
        currsize = axes2.OuterPosition;
        axes2.OuterPosition = [currsize(1:2) 950 420];
        axes2.Units = 'Normalized';

        % Title, axis labels and legend
        axes2.FontSize = 14;
        xlabel(axes2,[x_axis ' wavelength (cm^{-1})'],'FontSize',12,'FontWeight','bold');
        ylabel(axes2,'2D signal (a.u.)','FontSize',12,'FontWeight','bold')
        title(axes2,plot_title,'FontSize',10)
        axis(axes2,'tight');

        % Create legend
        legend(axes2,'show');
        legend(axes2,'boxoff')
        legend(axes2,'Location','eastoutside')

        % Add zero line
        hline = yline(axes2,0,'HandleVisibility','off');
        hline.Color = [0.5 0.5 0.5];
            
    case 'Integrate along probe axis'
        switch interactivemode
            case 'Off'
                % Ask user for pump range (default is same as probe range)
                defaults          = {num2str([ProbeAxis(1) ProbeAxis(end)])};
                SelTraces = inputdlg('Enter the pump range to show:',...
                             'Input pump wavenumber range to show:', [1 60],defaults);
                SelTraces = str2num(SelTraces{:});
            case 'On'
                SelTraces = [ProbeAxis(1) ProbeAxis(end)];
        end
        
        % Get the indices of the pump wavelengths
        pump_indexes    = findClosestId2Val(PumpAxis{1,1},SelTraces);
        % Initialize variables
        data            = zeros((max(pump_indexes)-min(pump_indexes)+1),Ndelays);
        % Get the data to be plotted
        for m=1:Ndelays
            data(:,m)   = sum(PROC_2D_DATA{m,1}(min(pump_indexes):max(pump_indexes),:),2);
        end
        
        % Create a new figure
        fh = figure();

        % Define the axes
        axes2 = axes(fh);

        % Plot the data
        cmap=colormap(axes2,othercolor('Mrainbow',Ndelays));
        for m=1:Ndelays
           plot(axes2,PumpAxis{1,1}(min(pump_indexes):max(pump_indexes)),data(:,m),'-','LineWidth',2,'MarkerSize',2,'color',cmap(m,:),'DisplayName',[num2str(t2delays(m),'%.3g') ' ps']);
           hold(axes2,'on');
        end

        % Save plot details
        plot_title  = '2D SIGNAL, INTEGRAL ACROSS ALL PROBE WAVELENGTHS';
        x_axis      = 'Pump';

        %%% Nice formatting
        % Size and colour
        fh.Color = [1 1 1];
        currsize = fh.Position;
        fh.Position = [currsize(1:2) 900 425];
        axes2.Units = 'Pixels';
        currsize = axes2.OuterPosition;
        axes2.OuterPosition = [currsize(1:2) 950 420];
        axes2.Units = 'Normalized';

        % Title, axis labels and legend
        axes2.FontSize = 14;
        xlabel(axes2,[x_axis ' wavelength (cm^{-1})'],'FontSize',12,'FontWeight','bold');
        ylabel(axes2,'2D signal (a.u.)','FontSize',12,'FontWeight','bold')
        title(axes2,plot_title,'FontSize',10)
        axis(axes2,'tight');

        % Create legend
        legend(axes2,'show');
        legend(axes2,'boxoff')
        legend(axes2,'Location','eastoutside')

        % Add zero line
        hline = yline(axes2,0,'HandleVisibility','off');
        hline.Color = [0.5 0.5 0.5];
        
    case 'Along several pump WL'
        switch interactivemode
            case 'On'
                % Select points interactively, to finish hit RETURN
                dataStruct = SelectTraces(dataStruct,0); 
                % Separate the values into pump and probe, according to how the data is plotted    
                switch plot_pumpdirection
                    case 'Horizontal'
                        pump_search    = dataStruct.SelTraces(:,1); 
                    case 'Vertical'
                        pump_search    = dataStruct.SelTraces(:,2); 
                end
            case 'Off'
                dataStruct.SelTraces = inputdlg('Enter the pump wavenumbers to plot:',...
                     'Input pump wavenumbers to plot:', [1 60]);
                pump_search = str2num(dataStruct.SelTraces{:})';
        end
        
        % Get the data to be plotted
        % Get the indices of the probe wavelengths
        pump_indexes        = unique(findClosestId2Val(PumpAxis{PopDelay,1},pump_search));
        probe_indexes       = 1:length(ProbeAxis);
        L                   = length(pump_indexes);
        % Initialise variables
        data            = zeros(length(ProbeAxis),L);
        % Get the data
        for p=1:L
            data(:,p)       = transpose(PROC_2D_DATA{PopDelay,1}(pump_indexes(p),probe_indexes));
        end
        
        if normalise
            data            = data./max(abs(data),[],1);
        end
        
        % Create a new figure
        fh = figure();

        % Define the axes
        axes2 = axes(fh);

        % Plot the data
        cmap=colormap(axes2,othercolor('Mrainbow',L));
        for p=1:L
           plot(axes2,ProbeAxis,data(:,p),'-','LineWidth',2,'MarkerSize',2,'color',cmap(p,:),'DisplayName',[num2str(PumpAxis{PopDelay,1}(pump_indexes(p)),'%.5g') ' cm^{-1}']);
           hold(axes2,'on');
        end

        % Save plot details
        plot_title  = ['CUTS ACROSS PUMP WAVELENGTHS AT t_2 = ' num2str(t2delays(PopDelay),'%.3g') ' ps'];
        x_axis      = 'Probe';

        %%% Nice formatting
        % Size and colour
        fh.Color = [1 1 1];
        currsize = fh.Position;
        fh.Position = [currsize(1:2) 900 425];
        axes2.Units = 'Pixels';
        currsize = axes2.OuterPosition;
        axes2.OuterPosition = [currsize(1:2) 950 420];
        axes2.Units = 'Normalized';

        % Title, axis labels and legend
        axes2.FontSize = 14;
        xlabel(axes2,[x_axis ' wavelength (cm^{-1})'],'FontSize',14,'FontWeight','bold');
        ylabel(axes2,'2D signal (a.u.)','FontSize',14,'FontWeight','bold')
        title(axes2,plot_title,'FontSize',10)
        axis(axes2,'tight');

        % Create legend
        legend(axes2,'show');
        legend(axes2,'boxoff')
        legend(axes2,'Location','eastoutside')

        % Add zero line
        hline = yline(axes2,0,'HandleVisibility','off');
        hline.Color = [0.5 0.5 0.5];
end
