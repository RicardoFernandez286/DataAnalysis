function plot3D(handles,Zdata,where)
% Get everything from handles structure
datafilename = handles.datafilename;
delays = handles.delays;
cmprobe = handles.cmprobe;
timescale = handles.timescale;
rawcorr = handles.rawcorr;
plotranges = handles.plotranges;
% Parse the plot ranges
delaylim = plotranges(1:2); %delaylim = [Min Max]
WLlim = plotranges(3:4); %WLlim = [Min Max]
Zminmax = plotranges(5); %Zminmax = max(|signal|) -> symmetric range
minabs = plotranges(7);
maxabs = plotranges(8);

% Do the plot
axes(where)
surf(cmprobe, delays, Zdata,'LineStyle','-','EdgeColor',[0.5 0.5 0.5],'EdgeAlpha','0.25','FaceColor','interp')
xlim(WLlim); ylim(delaylim);
colorbar;
% Make the plot format nice

% ORIGINAL COLOR DESIGN
% colormap(othercolor('BuDRd_18'));
% caxis([-Zminmax,Zminmax]);

% NEW COLOR DESIGN
% Define how many colours
n_tot = 500;
n_whites = floor(handles.percentwhites*n_tot/100);
percent_red = abs(maxabs/(abs(maxabs-minabs)));
percent_blue = 1 - percent_red;
n_blues = n_tot*percent_blue;
n_reds = n_tot*percent_red;

% Load the colormaps
map_red = othercolor('YlOrRd9',n_reds);
%map_blue = flipud(othercolor('PuBu9',n_blues));
ylord = othercolor('YlOrRd9',n_blues);
map_blue = flipud([ylord(:,3) ylord(:,2) ylord(:,1)]);
map_white = ones(n_whites,3);

% Renormalize the colour scale to make the white more striking
for j=1:3
    map_blue(:,j) = map_blue(:,j)/(max(map_blue(:,j)));
    map_red(:,j) = map_red(:,j)/(max(map_red(:,j)));
end

% Join the colormap "parts"
map = [map_blue;map_white;map_red];
colormap(map);
shading flat;
caxis([minabs,maxabs])

% Show the colorbar
hcb=colorbar;
title(hcb,{'\DeltaAbs';'(mOD)'});
set(gca,'FontSize',12);
title({datafilename;[rawcorr,' DATA';'']},'Interpreter','none','FontWeight','bold')
xlabel(['Wavenumbers',' (cm^{-1})'],'FontSize',13,'FontWeight','bold');
ylabel(['Delays',' (',timescale,')'],'FontSize',13,'FontWeight','bold');
zlabel('\DeltaAbs (mOD)');
view(-150,20)
hline = refline(0,0); hline.Color = 'black';
set(where,'yscale',handles.linlog)