m=4;
n=4;

N_phases    = m*n;
phase       = linspace(-pi,pi,N_phases);

fh = figure;
fh.Color = [1 1 1];
fh.Units = 'normalized';
fh.OuterPosition = [0 0 1 1];

for i=1:N_phases
    ax = subplot(m,n,i);
    Test_lab3_phasing(ax,phase(i))
    title(ax,['Phase: ' num2str(phase(i)) ' rad'],'FontSize',10);
end