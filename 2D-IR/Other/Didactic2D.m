fh = figure(1);
clf(fh)
ax1 = axes('parent',fh);

% Plot parameters
MinWL       = 1880;
MaxWL       = 2100;
Npoints     = 128;

Ncontours   = 40;
n_whites    = 2;

ContourLineColor = 0.65*[1 1 1];
ContourLineStyle = '-';
LineWidth   = 0.5;

LegendStyle = 'Both'; % 'PumpProbe' 'Omega' 'Both'

show1Dtop   = 1;
showPPside  = 0;
showContours= 1;

%% Calculate and plot

% Parameters
x1  = 1940;     % Centre of Peak 1
x2  = 2050;     % Centre of Peak 2

a1  = 20;       % Anharmonicity of Peak 1
a2  = 40;       % Anharmonicity of Peak 2

s1  = 11;       % Sigma of Peak 1
s2  = 12;       % Sigma of Peak 2

c1  = 0.4;      % Tilt of Peak 1 (0 <= c < 1)
c2  = 0.2;      % Tilt of Peak 2 (0 <= c < 1)
c12 = 0.1;      % Tilt of Cross-peak 1->2 (0 <= c < 1)
c21 = 0.1;      % Tilt of Cross-peak 2->1 (0 <= c < 1)

GSB1 = -1;      % Amplitude of GSB Peak 1
ESA1 = +1;      % Amplitude of ESA/SE Peak 1

GSB2 = -0.35;   % Amplitude of GSB Peak 2
ESA2 = +0.35;   % Amplitude of ESA/SE Peak 2

XPK12= 1;       % Relative amplitude of Xpeak 1->2
XPK21= 1;       % Relative amplitude of Xpeak 2->1

% Calculate preliminary things
x       = linspace(MinWL,MaxWL,Npoints);
y       = x;
[X,Y]   = meshgrid(x,y);

if show1Dtop == 1
    ax2 = axes('parent',fh);
end

if showPPside == 1
    ax3 = axes('parent',fh);
end

% Calculate the peaks
ZGSB    = G2Dc(X,Y,x1,x1,s1,s1,c1,GSB1/sqrt(1-c1))+G2Dc(X,Y,x2,x2,s2,s2,c2,GSB2/sqrt(1-c2));
ZESA    = G2Dc(X,Y,x1,x1-a1,s1,s1,c1,ESA1/sqrt(1-c1))+G2Dc(X,Y,x2,x2-a2,s2,s2,c2,ESA2/sqrt(1-c2));
ZXP12   = G2Dc(X,Y,x1,x2,s1,s2,c12,-XPK12*GSB1*GSB2/sqrt(1-c12))+G2Dc(X,Y,x1,x2-a2,s1,s2,c12,XPK12*ESA1*ESA2/sqrt(1-c12));
ZXP21   = G2Dc(X,Y,x2,x1,s2,s1,c21,-XPK21*GSB1*GSB2/sqrt(1-c21))+G2Dc(X,Y,x2,x1-a1,s2,s1,c21,XPK21*ESA1*ESA2/sqrt(1-c21));  

Z = ZGSB + ZESA + ZXP12 + ZXP21;

% Plotting stuff
switch LegendStyle
    case 'PumpProbe'
        PP_x    = 'Pump frequency (cm^{-1})';
        PP_y    = 'Probe frequency (cm^{-1})';
    case 'Omega'
        PP_x    = '\omega_{1} (cm^{-1})';
        PP_y    = '\omega_{3} (cm^{-1})';
    case 'Both'
        PP_x    = 'Pump frequency = \omega_{1} (cm^{-1})';
        PP_y    = 'Probe frequency = \omega_{3} (cm^{-1})';
end

max_cut         = max(abs(Z(:)));
min_cut         = -max_cut;
plot_contours   = linspace(min_cut,max_cut,Ncontours+1);
plot_outlines   = linspace(min_cut,max_cut,Ncontours/2);

%% Plot
hold(ax1,'on')
contourf(ax1,X,Y,Z,plot_contours,'LineColor','flat','LineStyle',ContourLineStyle,'LineWidth',LineWidth);
if showContours == 1
    contour(ax1,X,Y,Z,plot_outlines,'LineColor',ContourLineColor,'LineStyle',ContourLineStyle,'LineWidth',LineWidth*2);
end
hold(ax1,'off')

cmap    = darkb2r(min_cut,max_cut,Ncontours,n_whites);
colormap(ax1,cmap);
caxis(ax1,[min_cut,max_cut]);

diagline = refline(ax1,1,0);
diagline.Color = [0 0 0];
diagline.LineWidth = 0.5;

xlabel(ax1,PP_x,'FontWeight','bold','FontSize',14);
ylabel(ax1,PP_y,'FontWeight','bold','FontSize',14);

ax1.PlotBoxAspectRatio     = [1 1 1];
ax1.DataAspectRatio        = [1 1 1];
ax1.DataAspectRatioMode    = 'manual';
ax1.PlotBoxAspectRatioMode = 'manual';
ax1.Layer                  = 'top';
ax1.FontSize               = 14;
box(ax1,'on');

%% Beautify
if show1Dtop == 1
    TotalSum = -sum(ZGSB,2)./max(abs(sum(ZGSB,2)));
    plot(ax2,x,TotalSum,'Color','k','LineWidth',2);
    yl2 = yline(ax2,0);
    yl2.Color = [0.5 0.5 0.5];
    
    axis(ax2,'tight');
    ylim(ax2,[-0.1 1.1]);
    
    ylabel(ax2,'Abs (a.u.)','FontWeight','bold','FontSize',14);
    ax2.XTickLabel = [];

    ax2.Layer                  = 'top';
    ax2.FontSize               = 14;
    box(ax2,'on');
    
    ax1.Position    = [0.18    0.03    0.75    0.75];
    ax2.Position    = [0.18    0.72    0.75    0.25];
    fh.Position     = [1200 246 520 650];
    ax1.TickLength  = [0.02 0.01];
    ax2.TickLength  = [0.02 0.01];
else
    fh.Position     = [680 246 520 520];
    ax1.Position    = [0.18    0.18    0.75    0.75];
    ax1.TickLength  = [0.02 0.01];
end