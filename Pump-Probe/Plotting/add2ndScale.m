fh = gcf;
ax1 = gca;

if fh.Position(3) == 700
    fh.Position(3)  = 800;
    fh.Position(4)  = 425;
    fh.Color = 'w';
    title(ax1,'');
end

ax1.Units = 'pixels';
ax1.Position = [80 75 675 290];
ax1.Units = 'normalized';

FontSize = ax1.FontSize;
L = length(ax1.Children);

nm_axis = ax1.Children(1).XData;
cm_axis = 1e4./nm_axis;

for n=1:L
    ax1.Children(n).XData =  cm_axis;
end

ax1.XDir = 'reverse';

XL_cm = xlim(ax1);
XL_nm = 1e4./XL_cm;

if length(fh.Children) < 3
    ax_top = axes('parent',fh);
else
    ax_top = fh.Children(3);
end

cla(ax_top);

xlabel(ax1,'Wavenumbers (10^{3} cm^{-1})','FontSize',FontSize,'FontWeight','bold')
xlabel(ax_top,'Wavelength (nm)','FontSize',FontSize-2,'FontWeight','bold')
xlim(ax_top,[min(XL_nm), max(XL_nm)]);

ax_top.Color = 'none';
ax_top.Position = ax1.Position;
ax_top.YAxisLocation = 'right';
ax_top.XAxisLocation = 'top';
ax_top.Box = 'off';
ax_top.Color = 'white';
ax_top.FontSize  = FontSize-2;

linkaxes([ax1,ax_top],'y');
ax_top.YTickLabels = [];

ax1.Box = 'off';
ax1.Color = 'none';
axes(ax1);

XT_L = ax_top.XTickLabel;
XT_N = ax_top.XTick;

xlim(ax_top,[min(XL_cm), max(XL_cm)]);
ax_top.XTick = sort(1e4./XT_N);
ax_top.XTickLabel = flip(XT_L);
ax_top.XDir = 'reverse';

linkaxes([ax1,ax_top],'x');

ax_top.TickLength = ax1.TickLength;