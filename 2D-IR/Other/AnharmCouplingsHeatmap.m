% function AnharmCouplingsHeatmap

Nmodes = 4;
ModeLabels = ["A'(2)" "A''" "A'(1)" "C{\equiv}C"];

molec = 'NMe2';
cmapID = 'Blues5';
NanCol  = 1*[1 1 1];
LineCol = 1*[1 1 1];
LineWidth = 1.5;
LineStyle = 'solid'; % 'solid' | 'dashed' | 'dotted' | 'dashdot' |'none'

PlotW1W3  = 1; % 0 = energy ordering, 1 = 2D-IR like
plotTriang= -1; % 0 = plot everything; 1 = upper; -1 = lower;
hideDiag  = 0;

% freqs = [...
% 2002.62
% 2102.86
% 2226.10
% 3986.05
% 4151.60
% 4365.49
% 4078.34
% 4228.40
% 4327.87
% ];

%% Frequency input
switch molec
    case 'NMe2'
freqs = [...
1919.52
1942.86
2047.64
2214.33
3825.60
3874.80
4087.11
4410.94
3855.66
3953.72
4133.74
3971.61
4157.09
4261.84
];

    case 'Ph'
freqs = [...
1922.71
1945.21
2049.34
2254.69
3831.83
3879.38
4090.51
4491.78
3861.38
3958.33
4177.40
3975.65
4199.89
4303.99
];
end

%% Sort the data into the one and two exciton blocks
NH2 = Nmodes + nchoosek(Nmodes,2);
ExpectedN = Nmodes + NH2;

if length(freqs) ~= ExpectedN
    disp('Wrong number of modes')
    return
end

fund    = freqs(1:Nmodes);
ovcom   = freqs(Nmodes+1:end);

[ii,jj] = ndgrid(1:Nmodes);

H1 = fund*ones(Nmodes,1)'+(fund*ones(Nmodes,1)')';
H2 = zeros(Nmodes);

H2 = diag(ovcom(1:Nmodes),0);
H2(ii>jj) = ovcom(Nmodes+1:end);
H2 = triu(H2.',1) + tril(H2);

anh = H1-H2;

%% Plot and beautify
fh = figure(1);
fh.Color='w';
clf(fh);
fh.Position(3:4) = [509 420];

switch plotTriang % Plot a triangular matrix
    case -1
        if PlotW1W3 == 1
            anh(ii>jj) = nan; 
        else
            anh(ii<jj) = nan; 
        end
    case 1
        if PlotW1W3 == 1
            anh(ii<jj) = nan; 
        else
            anh(ii>jj) = nan; 
        end
end

if hideDiag == 1
    anh(ii==jj) = nan; % To remove diagonal
end

if PlotW1W3 == 1
    anh = flipud(anh);
end

h_map = heatmap(fh,anh,'MissingDataLabel','','CellLabelFormat','%0.1f','MissingDataColor',NanCol,'FontSize',16);
colormap(othercolor(cmapID,Nmodes*5));

warning('off','MATLAB:structOnObject');

hHeatmap = struct(h_map).Heatmap;
hGrid = struct(hHeatmap).Grid;
hGrid.ColorData = uint8(256*[LineCol 0.5]');
hGrid.LineWidth = LineWidth;
hGrid.LineStyle = LineStyle;

h_map.XData = ModeLabels;
if PlotW1W3 == 1
    h_map.YData = fliplr(ModeLabels);
else
    h_map.YData = ModeLabels;
end

%% Make it square

% Temporarily change axis units 
originalUnits = h_map.Units;  % save original units (probaly normalized)
h_map.Units = 'centimeters';  % any unit that will result in squares
% Get number of rows & columns
sz = size(h_map.ColorData); 
% Change axis size & position;
originalPos = h_map.Position; 
% make axes square (not the table cells, just the axes)
h_map.Position(3:4) = min(h_map.Position(3:4))*[1,1]; 
if sz(1)>sz(2)
    % make the axis size more narrow and re-center
    h_map.Position(3) = h_map.Position(3)*(sz(2)/sz(1)); 
    h_map.Position(1) = (originalPos(1)+originalPos(3)/2)-(h_map.Position(3)/2);
else
    % make the axis size shorter and re-center
    h_map.Position(4) = h_map.Position(4)*(sz(1)/sz(2));
    h_map.Position(2) = (originalPos(2)+originalPos(4)/2)-(h_map.Position(4)/2);
end
% Return axis to original units
h_map.Units = originalUnits; 

% % Temporarily change figure units 
% originalUnits_fig = fh.Units;  % save original units (probaly normalized)
% originalUnits_ax = h_map.Units; 
% fh.Units = 'centimeters';  % any unit that will result in squares
% h_map.Units = 'Normalized'; 
% % Get number of rows & columns
% sz = size(h_map.ColorData); 
% % make axes square (not the table cells, just the axes)
% h_map.Position(3:4) = min(h_map.Position(3:4))*[1,1]; 
% fh.InnerPosition(3:4) = min(fh.InnerPosition(3:4))*[1,1]; 
% % Change figure position;
% if sz(1)>sz(2)
%     % make the figure size more narrow
%     fh.InnerPosition(3) = fh.InnerPosition(3)*(sz(2)/sz(1)); 
% else
%     % make the figure size shorter
%     fh.InnerPosition(4) = fh.InnerPosition(4)*(sz(1)/sz(2)); 
% end
% % return original figure units
% fh.Units = originalUnits_fig; 
% h_map.Units = originalUnits_ax; 

% % Put an axis over top of the heatmap
% % The axis will not be visible but will cover all area 
% % to the right of the heatmap. 
% hmp = h_map.Position; 
% cbax = axes('Position',[sum(hmp([1,3])), 0, 1-sum(hmp([1,3])), 1],...
%     'XTick',[], 'YTick', [], 'Color',fh.Color);
% cbax.XAxis.Visible = 'off';
% cbax.YAxis.Visible = 'off';
% % Set the new axis color map to match the 
% % heatmap's colormap
% cbax.Colormap = h_map.Colormap; 
% % Add a colorbar the same vertical position as the heatmap
% cbh = colorbar(cbax,'Position',[.90, hmp(2), .03, hmp(4)]); 
% % Set the limits to 0:1 and set the ticks so that they 
% % are in the center of each bar
% cbh.Limits = [0,1]; 
% nColors = size(h_map.Colormap,1); 
% cbh.Ticks = (1/nColors : 1/nColors : 1) - 1/nColors/2; 
% % Set the new labels for each tick
% cbh.TickLabels = 0:nColors-1; 
% % Set the colorbar fontsize to the same value as heatmap fontsize
% cbh.FontSize = h_map.FontSize;

% % Temporarily change figure units 
% originalUnits_fig = fh.Units;  % save original units (probaly normalized)
% originalUnits_ax = h_map.Units; 
% fh.Units = 'centimeters';  % any unit that will result in squares
% h_map.Units = 'Normalized'; 
% % Get number of rows & columns
% sz = size(h_map.ColorData); 
% % make axes square (not the table cells, just the axes)
% h_map.Position(3:4) = min(h_map.Position(3:4))*[1,1]; 
% fh.InnerPosition(3:4) = min(fh.InnerPosition(3:4))*[1,1]; 
% % Change figure position;
% if sz(1)>sz(2)
%     % make the figure size more narrow
%     fh.InnerPosition(3) = fh.InnerPosition(3)*(sz(2)/sz(1)); 
% else
%     % make the figure size shorter
%     fh.InnerPosition(4) = fh.InnerPosition(4)*(sz(1)/sz(2)); 
% end
% % return original figure units
% fh.Units = originalUnits_fig; 
% h_map.Units = originalUnits_ax;


%% Plot a simulated 2D-IR spectrum with arbitrary lineshapes and intensities
MinWL   = (min(H1(:)))/2-40;
MaxWL   = (max(H1(:)))/2+30;
Npoints = 256;

x       = linspace(MinWL,MaxWL,Npoints);
y       = x;
[X,Y]   = meshgrid(x,y);

Z = zeros(size(X));

sx   = 8;
sy   = sx;
c1   = 0;
int  = [1 1 1 1];

for i=1:Nmodes
    for j=1:Nmodes
        x1 = fund(i);
        x3 = fund(j);
        Anh = H1-H2;
        a1 = Anh(i,j);
        ESA1 = (int(i)*int(j)).^2;
        GSB1 = -ESA1;
        % Calculate the peaks
        ZGSB    = G2Dc(X,Y,x1,x3,sx,sy,c1,GSB1/sqrt(1-c1));
        ZESA    = G2Dc(X,Y,x1,x3-a1,sx,sy,c1,ESA1/sqrt(1-c1));

        Z = Z + (ZGSB + ZESA)';
    end
end

fh = figure(2);
fh.Color = 'w';
clf(fh)
ax1 = axes('parent',fh);

PP_x    = '\omega_{1} (cm^{-1})';
PP_y    = '\omega_{3} (cm^{-1})';
percentScale = 100;
Ncontours = 40;
plot_skip = 4;
showContours = 1;

max_cut         = percentScale/100*max(abs(Z(:)));
min_cut         = -max_cut;
plot_contours   = linspace(min_cut,max_cut,Ncontours+1);
plot_contours   = plot_contours(1:end-1);

% negCont         = decim(linspace(min_cut,0,round((Ncontours+1)/2)),plot_skip,'mean');
% posCont         = decim(linspace(0,max_cut,round((Ncontours+1)/2)),plot_skip,'mean');
% plot_outlines   = [negCont; posCont];
plot_outlines   = linspace(min_cut,max_cut,round(Ncontours/plot_skip));

% Plot
cla(ax1);

contourf(ax1,X,Y,Z',plot_contours,'LineColor','none','LineStyle','-','LineWidth',0.1);
hold(ax1,'on')
if showContours == 1
    contour(ax1,X,Y,Z',plot_outlines,'LineColor',[0 0 0],'LineStyle','-','LineWidth',0.1);
end
hold(ax1,'off')

cmap    = darkb2r(min_cut,max_cut,Ncontours,2);
colormap(ax1,cmap);
caxis(ax1,[min_cut,max_cut]);

diagline = refline(ax1,1,0);
diagline.Color = [0 0 0];
diagline.LineWidth = 1;

xlabel(ax1,PP_x,'FontWeight','bold','FontSize',14);
ylabel(ax1,PP_y,'FontWeight','bold','FontSize',14);

ax1.PlotBoxAspectRatio     = [1 1 1];
ax1.DataAspectRatio        = [1 1 1];
ax1.DataAspectRatioMode    = 'manual';
ax1.PlotBoxAspectRatioMode = 'manual';
ax1.Layer                  = 'top';
ax1.FontSize               = 18;
box(ax1,'on');
fh.Position(3:4) = [643 551];
colorbar;