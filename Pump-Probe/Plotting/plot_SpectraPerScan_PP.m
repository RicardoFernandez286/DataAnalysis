function plot_SpectraPerScan_PP(app,dataStruct)

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
            dataStruct.SelTraces = yp;
        end
        hold(newaxes,'off')
        
        % Delete the extra figure
        close(fh)
    case 'Off'
        dataStruct.SelTraces = inputdlg('Enter ONE delay to plot:',...
             'Input desired delay:', [1 20]);
        dataStruct.SelTraces = str2double(char(dataStruct.SelTraces));
        if isempty(dataStruct.SelTraces)
            return
        end
end
% Once we have the point, look for the index
% Get the values to find in the TIME vector
value = dataStruct.SelTraces; % column 1= Wavelength = chosen x coordinate
% k = index of the closest match for the selected trace
k = findClosestId2Val(dataStruct.delays,value);
switch app.PP_DataFormat.Value
    case 'Pump-Probe [Lab 2]' % Lab 2
        % The first scan has no _n.csv termination
            Scan1 = char(strcat(tempdir,filesep,datafilename,'_signal.csv'));
            TempScanData = csvread(Scan1);
        if Nscans > 1
            ScanData(:,1) = TempScanData(k,:);
            % Load scans 2 to Nscans
            for i=2:Nscans
                ScanNr = strcat(tempdir,filesep,datafilename,'_signal_',num2str(i-1),'.csv');
                ScanNrch = char(ScanNr);
                TempScanData = csvread(ScanNrch);
                ScanData(:,i) = TempScanData(k,:);
                caption(i) = string(['Scan ' num2str(i)]);   
            end
            caption(1)= "Scan 1";
            clear TempScanData;
        else
            caption = "Scan 1";
            ScanData(:,1) = TempScanData(:,k);
        end
        ScanData = transpose(ScanData);
    case {'Pump-Probe [Lab 1 & Lab 4]','Pump-Probe [TRUVIS]'} % Lab 1 & Lab 4
        % The first scan has no _n.csv termination
        ending          = '_sp0_sm0_du0_';
        Nscans          = dataStruct.Nscans;
        RawDelays       = csvread([app.rootdir filesep datafilename filesep datafilename '_delays.csv']);
        ScanData        = zeros(Nscans,length(dataStruct.cmprobe{DET}));
        for i=1:Nscans
            ScanNr              = strcat(tempdir,filesep,datafilename,'_signal',ending,num2str(i-1),'.csv');
            ScanNrch            = char(ScanNr);
            TempScanData        = csvread(ScanNrch);
            [~,TempScanData]    = RemoveDuplicates(RawDelays,TempScanData);
            % Store it
            ScanData(i,:)       = TempScanData(k,:);
            % Caption it
            caption(i)          = string(['Scan ' num2str(i)]);   
        end
        clear TempScanData;
    case 'UniGE TA'
        ScanData = squeeze(dataStruct.scandata{1}(k,:,:))';
end

% Ask the user which scans to use if >15 scans
if Nscans >= 15
    plotscans   = inputdlg(['Please indicate which scans you want to plot (max: ' num2str(Nscans) '): '],'Select scans to plot');
    if strcmpi(plotscans,'all')
        plotscans  = 1:Nscans;
    else
        plotscans  = str2num(plotscans{:});
    end
    ScanData    = ScanData(plotscans,:);
    Nplots      = length(plotscans);
else
    Nplots      = Nscans;
    plotscans   = 1:Nscans;
end

% Bin the scans if BinScans > 1
if BinScans > 1
    i=1; s=1;
    BinScanData=zeros(ceil(Nplots/BinScans),length(dataStruct.cmprobe{DET}));
    while s <= Nplots
        snew = s+BinScans;
        if snew > Nplots % Last bin
            m=s:Nplots;
            BinCaption(i)=string([num2str(plotscans(s)) ' to ' num2str(plotscans(Nplots))]);
        else
            m=s:snew-1;
            BinCaption(i)=string([num2str(plotscans(s)) ' to ' num2str(plotscans(snew-1))]);
        end
        BinScanData(i,:)=mean(ScanData(m,:),1);
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
        for j=1:Nplots
            ScanData(j,:) = ScanData(j,:)./max(abs(ScanData(j,:)),[],'all');
        end
end

% Create a new figure with consistent format
fh 				= figure();
fh.Position(3)  = 880;
fh.Position(4)  = 425;
fh.Color        = [1 1 1];

% Define the axes
axes2 = axes('Parent',fh);
axes(axes2);

% Plot the data
cm=colormap(othercolor('Mrainbow',Nplots));
for n=1:Nplots
   pl(n) = plot(axes2,dataStruct.cmprobe{DET},ScanData(n,:),'color',cm(n,:));
   hold on
end
% Make the lines of the first and last scans thicker
pl(1).LineWidth     = 1.5;
pl(end).LineWidth   = 1.5;

% Nice formatting
axes2.FontSize 		= 12;
title(axes2,{dataStruct.datafilename;['TR. SPECTRA PER SCAN AT ' num2str(dataStruct.delays(k),'%.3g') ' ' dataStruct.timescale]},'Interpreter','none','FontSize',12)
xlabel(axes2,dataStruct.probeunits,'FontSize',13,'FontWeight','bold')
ylabel(axes2,label,'FontSize',13,'FontWeight','bold')
axis(axes2,'tight');
hline       = yline(axes2,0,'HandleVisibility','off'); 
hline.Color = [0.5 0.5 0.5];
legend(axes2,caption);
legend(axes2,'boxoff');
legend(axes2,'Location','bestoutside');
axes2.Units     = 'pixels';
axes2.Position  = [75 60 675 320];
axes2.Units     = 'normalized';