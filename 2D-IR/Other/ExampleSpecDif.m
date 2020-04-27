function ExampleSpecDif
% Nc=20;
% C = linspace(0.9,0,Nc);
Nc=1;
C = 0.;
F(Nc) = struct('cdata',[],'colormap',[]);
DoMovie = 0;

fh = figure(1);
fh.Color = 'w';
clf(fh)
ax1 = axes('parent',fh);

% Plot parameters
MinWL       = 1950;
MaxWL       = 2040;
Npoints     = 256;

Ncontours   = 40;
n_whites    = 2;

ContourLineColor = 0.5*[1 1 1];
ContourLineStyle = '-';
LineWidth   = 0.5;

LegendStyle = 'Omega'; % 'PumpProbe' 'Omega' 'Both'

show1Dtop   = 1;
showPPside  = 0;
showContours= 1;
randomColors= 0;
show1Dinhom = 1;
show2DverticalLines = 0;
%% Calculate and plot

% Parameters
x1  = 2000;     % Centre of Peak 1

a1  = 20;       % Anharmonicity of Peak 1

s1  = 10;       % Sigma of Peak 1

% c1  = 0.8;      % Tilt of Peak 1 (0 <= c < 1)


GSB1 = -1;      % Amplitude of GSB Peak 1
ESA1 = +1;      % Amplitude of ESA/SE Peak 1

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

for k=1:Nc
    c1 = C(k);
% Calculate the peaks
ZGSB    = G2Dc(X,Y,x1,x1,s1,s1,c1,GSB1/sqrt(1-c1));
ZESA    = G2Dc(X,Y,x1,x1-a1,s1,s1,c1,ESA1/sqrt(1-c1));

Z = ZGSB + ZESA;

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
contourf(ax1,X,Y,Z,plot_contours,'LineColor','flat','LineStyle',ContourLineStyle,'LineWidth',LineWidth);
hold(ax1,'on')
if showContours == 1
    contour(ax1,X,Y,Z,plot_outlines,'LineColor',ContourLineColor,'LineStyle',ContourLineStyle,'LineWidth',LineWidth*2);
end
hold(ax1,'off')

cmap    = darkb2r(min_cut,max_cut,Ncontours,n_whites);
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
ax1.FontSize               = 14;
box(ax1,'on');


%% Beautify
if show1Dtop == 1
    TotalLinear = -sum(ZGSB,2)./max(abs(sum(ZGSB,2)));
    
    if show1Dinhom == 1    
        Ncolors = 8;
        Nskip = 1;
%         sH = s1*(Ncolors/2);
        sH = s1/4;
        
        Wx = linspace(1975,2025,Ncolors);
        cmap = flipud(othercolor('Mrainbow',Ncolors));
%         cmap = flipud(othercolor('Spectral10',Ncolors));
%         cmap = flipud(jet(Ncolors));
        
        if randomColors == 1
            cmap = cmap(randperm(length(cmap)),:);
        end
        
        hold(ax2,'on');
        for i=1:Nskip:Ncolors
            plot(ax2,x,exp(-((x-Wx(i))/sH).^2).*TotalLinear','Color',cmap(i,:),'LineWidth',1.5);
        end
              
        if show2DverticalLines == 1
            for i=1:Ncolors
                xline(ax1,Wx(i),'Color',cmap(i,:),'LineWidth',1.5);
            end
        end
    end
    yl2 = yline(ax2,0,'Color',[0.5 0.5 0.5],'LineWidth',2);
    
    plot(ax2,x,TotalLinear,'Color','k','LineWidth',2);
    hold(ax2,'off');
    
    axis(ax2,'tight');
    ylim(ax2,[-0.1 1.1]);
    
    ylabel(ax2,'Abs (a.u.)','FontWeight','bold','FontSize',14);
    ax2.XTickLabel = [];

    ax2.Layer                  = 'top';
    ax2.FontSize               = 14;
    box(ax2,'on');
    
    linkaxes([ax1 ax2],'x');
    
    xlim(ax1,[1960 2040])
    ylim(ax1,[1950 2030])
    ax1.Position    = [0.18    0.04    0.75    0.75];

    ax2.Position    = [0.18    0.74    0.75    0.25];
    fh.Position     = [1500 300 520 680];
    ax1.TickLength  = [0.02 0.01];
    ax2.TickLength  = [0.02 0.01];
    ax1.LineWidth   = 1;
    ax2.LineWidth   = 1;
    if DoMovie == 1
        drawnow;
        F(k) = getframe(fh);
    end
else
    fh.Position     = [680 246 520 520];
    ax1.Position    = [0.18    0.18    0.75    0.75];
    ax1.TickLength  = [0.02 0.01];
end
end

if DoMovie == 1
    % Write to the GIF File
    movie2gif(F,'specdif.gif','DelayTime',0.2,'LoopCount',Inf);
end