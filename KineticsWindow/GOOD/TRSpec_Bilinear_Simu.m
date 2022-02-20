function [WL,t,D,S] = TRSpec_Bilinear_Simu(MakePlots,NoisePercent)
% Description:  This script simulates time-dependent spectra for a given set of components given their spectra and 
%               the kinetic parameters needed to build a concentration matrix based on coupled first-order reactions
%
% Tested and implemented in MATLAB R2021b
% v1.0

%% Build S matrix
% Define a Wavelength axis
WL  = linspace(300,800,500);
nWL = length(WL);

% Define the spectra of each component (assuming Gaussian lineshapes for simplicity)
% Each species will contain up to 3 Gaussians (Row=species, Column=Gaussian subcomponent
A   = [0.8  0      0
       0.75 0.45   0
       0.2  0.15   0.3
       0.2  0.7    0];

fwhm= [80   1      1
       80   60     1
       50   50     60
       80   80     0];

W0  = [380  0      0
       500  650    0
       430  610    680
       600  720    0];

S   = zeros(4,nWL);

for i=1:4
    for j=1:3
        S(i,:) = S(i,:) + A(i,j).*exp(-((WL-W0(i,j))./fwhm(i,j)).^2);
    end
end

%% Build C matrix
k1  = 0.5;
k2f = 0.3;
k2b = 0.4;
k3  = 0.1;

K   = [-k1  0       0         0
	   k1   -k2f    k2b       0
	   0    k2f     -k3-k2b   0
	   0    0       k3        0];

nt  = 250;
t   = [0 logspace(-2,2,nt-1)]';
C0  = [1 0 0 0]';

C = kineticsKmat_simu(t,K,C0,1*MakePlots,[]);
%% Build transient spectra
D = C*S;
D = D + NoisePercent.*max(abs(D(:)))./100.*randn(nt,nWL);

if MakePlots == 1
%% Make diverse plots
fh=figure(1);  fh.Name = 'Concentration Profiles';

% Make a contour plot of the data (with X and Y lines to show where the traces were taken)
fh      = figure(2); fh.Name = 'Contour Plot';
fh.Color='w';
clf(fh)
ax      = axes('parent',fh);

maxPlot = max(D(:));
minPlot = min(D(:));
NCtrs   = 40;
ctrs    = linspace(minPlot,maxPlot,NCtrs);

contourf(ax,WL,t,D,ctrs,'EdgeColor','none')

colormap(ax,turbo(NCtrs+1))
ax.Layer        = 'top';
ax.YScale       = 'log';
ax.TickLength   = [0.025 0.025];

cb = colorbar(ax);
ylabel(cb,'Absorbance','FontWeight','bold');

ylabel(ax,'Time','FontWeight','bold')
xlabel(ax,'Wavelength (nm)','FontWeight','bold')
ax.FontSize = 14;

% Plot the transient spectra (cuts across the time axis)
TraceT      = [0 0.1 0.25 0.5 0.75 1 2.5 5 10 25 50 100];
TraceIdxT   = sort(unique(dsearchn(t,TraceT')));
nTt         = length(TraceIdxT);

cmapT = turbo(nTt+1);

fh = figure(3); fh.Name = 'Transient Spectra';
fh.Color='w';
clf(fh)
ax=axes('parent',fh);

hold(ax,'on')
for j=1:nTt
plot(WL,D(TraceIdxT(j),:),'LineWidth',2,'Color',cmapT(j,:),'DisplayName',[num2str(t(TraceIdxT(j)),3)])
% drawnow;
% pause(0.025)
end
hold(ax,'off')

xlabel(ax,'Wavelength (nm)','FontWeight','bold');
ylabel(ax,'Absorbance','FontWeight','bold');

lh = legend(ax,'show','Location','bestoutside','box','off');
title(lh,'Time (s)')
ax.Box          = 'on';
ax.TickLength   = [0.025 0.025];
ax.FontSize     = 14;

% Plot selected kinetic traces
TraceWL = [380 450 500 600 700 750];
TraceIdx= sort(unique(dsearchn(WL',TraceWL')));
nTWL    = length(TraceIdx);

cmapL = turbo(nTWL+1);

fh = figure(4); fh.Name = 'Kinetic Traces';
fh.Color='w';
clf(fh)
ax=axes('parent',fh);
hold(ax,'on');
for i=1:nTWL
    plot(ax,t,D(:,TraceIdx(i)),'LineWidth',2,'Color',cmapL(i,:),'DisplayName',[num2str(round(WL(TraceIdx(i)))) ' nm'])
end
hold(ax,'off');

xlabel(ax,'Time','FontWeight','bold');
ylabel(ax,'Absorbance','FontWeight','bold');

legend(ax,'show','Location','best','box','off')
ax.Box          = 'on';
ax.XScale       = 'log';
ax.TickLength   = [0.025 0.025];
ax.FontSize     = 14;


end