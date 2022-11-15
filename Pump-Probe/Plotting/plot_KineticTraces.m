function [figNum,k,ydata,yFit,res] = plot_KineticTraces(DataStruct,DET,data,overlay,figNum,Normalise,varargin)

%% Hardcoded Settings
FontSize = 16; % 18 for papers (single column figure), 14 for presentations (can vary)
ResPlotStyle = 'o';
%% Decide whether we want residuals or not
if ~isempty(varargin)
    resplot = varargin{1};
    if resplot ~= 1
        overlay = 0;
    end
end

% Empty variables in case they are not overwritten
res     = [];
yFit    = [];
%% Get the data to be plotted
% Get the values to find in the WAVELENGTH vector
value = DataStruct.SelTraces(:,1); % column 1= Wavelength
% k = index of the closest match for each selected trace
k = unique(findClosestId2Val(DataStruct.cmprobe{DET},transpose(value)));
L = length(k);
i = 0;
while i < L
    % Get the actual vectors from the Z matrix (one by one)
    index(i+1) = k(1,i+1); %#ok<*AGROW> 
    switch Normalise
        case 0
            ydata(:,i+1) = data(:,index(i+1));
            label = '\DeltaAbs (mOD)';
            caption(i+1,1) = string([num2str(DataStruct.cmprobe{DET}(index(i+1))) ' ' DataStruct.Xunits]);
        case 1
            maxval = max(data(:,index(i+1)));
            minval = min(data(:,index(i+1)));
            if abs(maxval) >= abs(minval)
                ydata(:,i+1) = data(:,index(i+1))./maxval;
                caption(i+1,1) = string([num2str(DataStruct.cmprobe{DET}(index(i+1))) ' ' DataStruct.Xunits]);
            else
                ydata(:,i+1) = data(:,index(i+1))./minval;
                caption(i+1,1) = string([num2str(DataStruct.cmprobe{DET}(index(i+1))) ' ' DataStruct.Xunits '\times -1']);
            end
            label = 'Normalised \DeltaAbs';
    end
    i=i+1;
end

% Check that I have a fit to plot, otherwise don't do it
if overlay == 2
    if ~isfield(DataStruct,'FitData') || isempty(DataStruct.FitData{DET})
        overlay = 0;
    else
        FitData = DataStruct.FitData{DET};
        switch Normalise
            case 0
                yFit    = FitData(:,k);
            case 1
                maxval  = max(data(:,k),[],1);
                minval  = min(data(:,k),[],1);
                sgn     = sign(abs(maxval)-abs(minval));
                posT    = find(~isnan(FitData(:,k)),1);
%                 posT    = DataStruct.delays > 0;
                ydata   = sgn.*data(:,k)./max(abs(data(posT:end,k)));
                yFit    = sgn.*FitData(:,k)./max(abs(data(posT:end,k)));
        end
        res = ydata - yFit;
    end
end

%% Make figure and axes
if overlay == 1
    h =  findobj('type','figure');
    n = length(h);
    if n > 0
        figNum = n;
    else
        figNum = 1;    
    end
    fh      = figure(figNum);
    
    if isempty(fh.Children)
        axes2   = axes('Parent',fh);
    else
        axes2   = fh.Children(end);
    end
    
    hold(axes2,'on')
    axes(axes2);
elseif overlay == 2
    % Create a new figure
    fh = figure;
    fh.Position([3 4]) = [890 540];
    % % Define the axes
    axes2 = axes('Parent',fh);
    axes(axes2);
    resaxes = axes('Parent',fh);

    axes2.Visible = 'off';
    resaxes.Visible = 'off';
else
    % Create a new figure
    fh = figure;
    fh.Position([3 4]) = [920 425];
    % % Define the axes
    axes2 = axes('Parent',fh);
    axes(axes2);
    
    % % Plot in two axes 
    % figure;
    % ax1 = axes('OuterPosition',[0 0 0.3 1],'Units','Normalized');
    % ax2 = axes('OuterPosition',[0.1 0 0.8 1],'Units','Normalized','YTickLabel',[]);
end

movegui(fh,'onscreen');
fh.Color        = [1 1 1];
%% Plot the data
cmap=colormap(axes2,(othercolor('Mrainbow',L)));

if overlay == 2
   hold(axes2,'on');
   hold(resaxes,'on');
    for n=1:L
       plot(axes2,DataStruct.delays,ydata(:,n),'o','MarkerSize',5,'MarkerEdgeColor',cmap(n,:),'MarkerFaceColor','none','HandleVisibility','off');
       plot(axes2,DataStruct.delays,yFit(:,n),'-','LineWidth',2,'Color',cmap(n,:),'DisplayName',[num2str(DataStruct.cmprobe{DET}(index(n)),'%.1f') ' ' DataStruct.Xunits]);
       if Normalise == 1
           plot(resaxes,DataStruct.delays,100.*res(:,n),ResPlotStyle,'MarkerSize',5,'Color',cmap(n,:),'MarkerEdgeColor',cmap(n,:),'MarkerFaceColor','none','DisplayName',[num2str(DataStruct.cmprobe{DET}(index(n)),'%.1f') ' ' DataStruct.Xunits]);
       else
           plot(resaxes,DataStruct.delays,res(:,n),ResPlotStyle,'MarkerSize',5,'Color',cmap(n,:),'MarkerEdgeColor',cmap(n,:),'MarkerFaceColor','none','DisplayName',[num2str(DataStruct.cmprobe{DET}(index(n)),'%.1f') ' ' DataStruct.Xunits]);
       end
       xline(resaxes,0,'HandleVisibility','off');
       yline(resaxes,0,'HandleVisibility','off');
    end
    hold(axes2,'off');
    hold(resaxes,'off');
else
    for n=1:L
       plot(axes2,DataStruct.delays,ydata(:,n),'LineWidth',2,'Marker','o','MarkerSize',2,'color',cmap(n,:),'DisplayName',[num2str(DataStruct.cmprobe{DET}(index(n)),'%.1f') ' ' DataStruct.Xunits]);
       hold(axes2,'on');
    end
    hold(axes2,'off');
end
%% Nice formatting
title(axes2,{DataStruct.datafilename;[DataStruct.rawcorr,' DATA - KINETICS';'']},'Interpreter','none','FontSize',12);
xlabel(axes2,['Delays',' (',DataStruct.timescale,')'],'FontSize',13,'FontWeight','bold');
ylabel(axes2,label,'FontSize',13,'FontWeight','bold');
axis(axes2,'tight');
%             if L == 1
%                 legend(axes2,[string(strcat(num2str(round(DataStruct.cmprobe(index),0)),[' ' DataStruct.Xunits])),''],'FontSize',12)
%             else
%                 legend(axes2,caption,'FontSize',12)
%             end
legend(axes2,'show')
legend(axes2,'boxoff')
legend(axes2,'Location','bestoutside')
axes2.XScale    = DataStruct.linlog;
yline(axes2,0,'Color',[0.5 0.5 0.5],'HandleVisibility','off');
xline(axes2,0,'Color',[0.5 0.5 0.5],'HandleVisibility','off');
axes2.Units     = 'normalized';
axes2.Layer     = 'top';
axes2.Box       = 'on';
axes2.TickLength=[0.015 0.025];
axes2.FontSize  = FontSize;

if overlay == 1
   Nplots = length(fh.Children(end).Children);
   cmap = othercolor('Mrainbow',Nplots);
   for j=1:Nplots
       fh.Children(end).Children(j).Color = cmap(j,:);
   end
   if Nplots > 1
       title(axes2,'Overlay plot');
   end
elseif overlay == 2
   resaxes.XScale    = DataStruct.linlog;
   linkaxes([resaxes,axes2],'x');
   
   xlim(resaxes,'tight');
   
   axes2.Units     = 'pixels';
   axes2.Position  = [75 60 675 320];
   resaxes.Units= 'pixels';
   axes2.Units  = 'pixels';

   axes2.Position([1 3]) = [80 640];
   resaxes.Position([1 3]) = axes2.Position([1 3]);
   resaxes.Position([2 4]) = [80 140];
   axes2.Position([2 4])   = [230 250];

   axes2.XTickLabel  = [];
   resaxes.Layer     = 'top';
   resaxes.Box       = 'on';
   resaxes.TickLength= [0.015 0.025];

   resaxes.FontSize  = FontSize;
   xlabel(resaxes,['Delays',' (',DataStruct.timescale,')'],'FontSize',FontSize,'FontWeight','bold');
   if Normalise == 1
        ylabel(resaxes,'% Res.','FontSize',FontSize,'FontWeight','bold');
   else
        ylabel(resaxes,'Res. (mOD)','FontSize',FontSize,'FontWeight','bold');
   end
   xlabel(axes2,[]);
   
   resaxes.Units     = 'normalized';
   axes2.Units       = 'normalized';
   
   axes2.Visible     = 'on';
   resaxes.Visible   = 'on';
end