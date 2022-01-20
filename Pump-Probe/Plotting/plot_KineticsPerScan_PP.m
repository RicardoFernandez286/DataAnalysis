function plot_KineticsPerScan_PP(app,dataStruct)

Nscans          = dataStruct.Nscans;
datafilename    = dataStruct.datafilename;
tempdir         = [app.rootdir filesep datafilename filesep 'temp'];
tempdirOUT      = [app.rootdir filesep datafilename 'temp'];

DET             = app.PP_ProbeDetSwitch.Value;
interactive     = app.PP_InteractivemodeSwitch.Value;
normalise       = app.PP_NormaliseCheckBox.Value;
BinScans        = dataStruct.BinScans;

% Check if the temp dir is outside (Lab 1)
if exist(tempdir,'dir') == 0
    tempdir     = tempdirOUT;
end

dataStruct.SelTraces    = [];
switch interactive
    case 'On'
        % Patch because ginput doesn't work on UIFigures
        fh = figure;
        fh.Position = app.PumpProbeAxes.Position;
        fh.Position(1:2) = fh.Position(1:2) + app.DataAnalysisGUI_UIFigure.Position(1:2);
        newaxes = axes('Parent',fh);
        UpdatePlot_PP(app,newaxes);

        % Select ONE point to compare the KINETICS
        hold(newaxes,'on')
        [xp,yp,button] = ginput(1);
        if isequal(button,1)
            plot(newaxes,xp,yp,'X','Linewidth',1.5,'Color','k')
            dataStruct.SelTraces = xp;
        end
        hold(newaxes,'off')
        
        % Delete the extra figure
        close(fh)
    case 'Off'
        dataStruct.SelTraces = inputdlg('Enter ONE probe wavelength to plot:',...
             'Input desired probe frequency:', [1 50]);
        dataStruct.SelTraces = str2double(dataStruct.SelTraces{:})';
        if isempty(dataStruct.SelTraces)
            return
        end
end

% Once we have the point, look for the index
% Get the values to find in the WAVELENGTH vector
value   = dataStruct.SelTraces; % column 1= Wavelength = chosen x coordinate
% k = index of the closest match for the selected trace
k       = findClosestId2Val(dataStruct.cmprobe{DET},value);
% Load the scan data depending on Lab 1 vs Lab 2
switch app.PP_DataFormat.Value
   case 'Pump-Probe [Lab 2]'
        % The first scan has no _n.csv termination
            Scan1           = [tempdir filesep datafilename '_signal.csv'];
            TempScanData    = csvread(Scan1);
            RawDelays       = csvread([app.rootdir filesep datafilename filesep datafilename '_delays.csv']);
        if Nscans > 1
            ScanData(:,1) = TempScanData(:,k);
            % Load scans 2 to Nscans
            for i=2:Nscans
                ScanNr          = [tempdir filesep datafilename '_signal_' num2str(i-1) '.csv'];
                ScanNrch        = char(ScanNr);
                TempScanData    = csvread(ScanNrch);
                ScanData(:,i)   = TempScanData(:,k);
                caption(i)      = string(['Scan ' num2str(i)]);   
            end
            caption(1)= "Scan 1";
            clear TempScanData;
        else
            caption = "Scan 1";
            ScanData(:,1) = TempScanData(:,k);
        end
    case {'Pump-Probe [Lab 1 & Lab 4]','Pump-Probe [TRUVIS]'}
        % The first scan has _0.csv termination, then _(n-1).csv
        ending      = '_sp0_sm0_du0_';
        RawDelays   = csvread([app.rootdir filesep datafilename filesep datafilename '_delays.csv']);
        ScanData    = zeros(length(unique(RawDelays)),Nscans);
        for i=1:Nscans
            ScanNr              = [tempdir filesep datafilename '_signal' ending num2str(i-1) '.csv'];
            ScanNrch            = char(ScanNr);
            TempScanData        = csvread(ScanNrch);
            [~,TempScanData]    = RemoveDuplicates(RawDelays,TempScanData);
            % Store it
            ScanData(:,i) = TempScanData(:,k);   
        end
%         clear TempScanData;
    case 'UniGE TA'
        ScanData = squeeze(dataStruct.scandata{1}(:,k,:));
end
% Ask the user which scans to use if >15 scans
if Nscans >= 15
    plotscans   = inputdlg(['Please indicate which scans you want to plot (max: ' num2str(Nscans) '): '],'Select scans to plot');
    if strcmpi(plotscans,'all')
        plotscans  = 1:Nscans;
    else
        plotscans  = str2num(plotscans{:});
    end
    ScanData    = ScanData(:,plotscans);
    Nplots      = length(plotscans);
else
    Nplots      = Nscans;
    plotscans   = 1:Nscans;
end

% Bin the scans if BinScans > 1
if BinScans > 1
    i=1; s=1;
    BinScanData=zeros(length(dataStruct.delays),ceil(Nplots/BinScans));
    while s < Nplots
        snew = s+BinScans;
        if snew > Nplots % Last bin
            m=[s:Nplots];
            BinCaption(i)=string([num2str(plotscans(s)) ' to ' num2str(plotscans(Nplots))]);
        else
            m=[s:snew-1];
            BinCaption(i)=string([num2str(plotscans(s)) ' to ' num2str(plotscans(snew-1))]);
        end
        BinScanData(:,i)=mean(ScanData(:,m),2);
        s=snew; i=i+1;
    end
    RawScanData = ScanData;
    ScanData    = BinScanData;
    Nplots      = ceil(Nplots/BinScans);
    caption     = BinCaption;
else
    for i=1:Nplots
        % Caption it
        caption(i) = string(['Scan ' num2str(plotscans(i))]);
    end
end

% Normalise the data if NormKin = 1
switch normalise
    case 0
        label = '\DeltaAbs (mOD)';
    case 1
        label = 'Normalised \DeltaAbs';
        j=0;
        while j < Nplots
            j=j+1;
            ScanData(:,j) = ScanData(:,j)./max(max(abs(ScanData(:,j))));
        end
end

% Create a new figure
fh = figure();
fh.Position(3)  = 920;
fh.Position(4)  = 425;
fh.Color        = [1 1 1];

% Define the axes
axes2 = axes('Parent',fh);
axes(axes2);

% Plot the data
cm=colormap(axes2,othercolor('Mrainbow',Nplots));
for n=1:Nplots
   plot(axes2,dataStruct.delays,ScanData(:,n),'color',cm(n,:),'DisplayName',caption(n),'HandleVisibility','off');
   plot(axes2,NaN,NaN,'color',cm(n,:),'LineWidth',3)
   hold(axes2,'on')
end
% Make the lines of the first and last scans thicker
axes2.Children(end).LineWidth   = 2*axes2.Children(end).LineWidth;
axes2.Children(1).LineWidth     = 2*axes2.Children(1).LineWidth;

% Nice formatting
data_title  = regexprep(dataStruct.datafilename, '\_', '\\\_');
title(axes2,{data_title;['KINETICS PER SCAN AT ' num2str(round(dataStruct.cmprobe{DET}(k),0)) ' ' dataStruct.Xunits]},'Interpreter','tex')
xlabel(axes2,['Delays',' (',dataStruct.timescale,')'],'FontSize',13,'FontWeight','bold')
ylabel(axes2,label,'FontSize',13,'FontWeight','bold')
axis(axes2,'tight');
legend(axes2,caption);
legend(axes2,'Location','bestoutside');
legend(axes2,'boxoff');
axes2.XScale = dataStruct.linlog;
hline = yline(axes2,0,'HandleVisibility','off');
hline.Color = [0.5 0.5 0.5];
axes2.FontSize  = 12;
axes2.Units     = 'pixels';
axes2.Position  = [75 60 675 320];
axes2.Units     = 'normalized';