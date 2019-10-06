folder = '\\idnetapp-chem.uzh.ch\g_chem_hamm$\Group\Ricardo\ParamOpt\T2 VAR';

folderlist = dir(folder);
foldernames = {folderlist.name}';
foldernames = foldernames(3:end);

Nplots = length(foldernames);

for i=1:Nplots
    z_str = strsplit(foldernames{i},'_');
    z(i) = str2double(z_str{2});
end

[~,idx] = sort(z);

fh= figure(1);
clf(fh);
ax=axes('parent',fh);
hold(ax,'on')


wb = waitbar(0,'Loading...');
Nstop = Nplots;
cmap = jet(Nstop);
for i=1:Nstop
    waitbar(i/Nplots,wb,['Loading ' num2str(i) ' of ' num2str(Nplots) ' ...']);
    spec = dlmread([folder filesep foldernames{idx(i)} filesep 'spec_lin.dat']);
    x = spec(:,1);
    y(:,i) = spec(:,2);
    ynorm(:,i) = spec(:,2)/max(spec(:,2));
    z_str = strsplit(foldernames{idx(i)},'_');
    plot(x,ynorm(:,i),'Color',cmap(i,:),'DisplayName',z_str{2})
%     drawnow;
end
legend(ax,'show');

plot(exp_spec(:,1)-2,exp_spec(:,3),'LineWidth',3,'Color','k');
xlim([1925,2075]);
delete(wb);

