%% Plot_nmInverse
nLines      = length(ax.Children);

Ticks_nm    = ax.XTick;

% Filter the ones in the red edge
cutoff = 451;
Ticks_nm    = [Ticks_nm(Ticks_nm<cutoff) min(round(Ticks_nm(Ticks_nm>=cutoff)/100)*100):100:max(Ticks_nm)];
Ticks_cm    = 1e4./Ticks_nm;

nTicks      = length(Ticks_nm);

% Labels_nm   = ax.XTickLabel;
% Labels_cm   = cell(size(Labels_nm));

Labels_nm   = cell(nTicks,1);
Labels_cm   = cell(nTicks,1);

for i=1:nTicks
    Labels_nm{i} = num2str(Ticks_nm(i));
    Labels_cm{i} = num2str(Ticks_cm(i));
end

XL_nm = xlim(ax);
XL_cm = sort(1e4./XL_nm);

for i=1:nLines
    ax.Children(i).XData = 1e4./ax.Children(i).XData;
end

ax.XTick        = flip(Ticks_cm);
ax.XTickLabel   = flip(Labels_nm);

xlim(ax,XL_cm);
ax.XDir = 'reverse';

ax.XTickLabelRotation = 0;