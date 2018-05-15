normalized = 1;
TitleOnOff = 'Off';

name400='ReNCS_ZrO2_100mM_QH2_dB_ACN_d3_2uJ_1029';
% name400='RePOH_sC-ZrO2-164.1_2µJ_MeCN_200mM_QH2_1546';
name0='ReNCS_ZrO2_dB_ACN_d3_2uJ_091348';

sig400=csvread([name400 filesep name400 '_signal.csv']);
sig0=csvread([name0 filesep name0 '_signal.csv']);

sig400 = sig400 - mean(sig400(1:5,:),1);
sig0 = sig0 - mean(sig0(1:5,:),1);

sig400 = sig400 - mean(sig400(:,:),2);
sig0 = sig0 - mean(sig0(:,:),2);

if normalized == 1
    diff = sig400./max(max(abs(sig400))) - sig0./max(max(abs(sig0)));
else
    diff = sig400 - sig0;
end

delays=csvread([name0 filesep name0 '_delays.csv']);
wavenumbers=csvread([name0 filesep name0 '_wavenumbers.csv']);


contourf(wavenumbers,delays,diff,30,'LineColor','flat')
Zminmax=max(max(abs(diff)));

% Original colormap style:
colormap(othercolor('BuDRd_18',50));
caxis([-Zminmax,Zminmax])
pepita=gca;
set(pepita,'yscale','log')

% Show the colorbar
hcb=colorbar;
hcb.FontSize=12;
hcb.LineWidth = 0.1;
hcb.TickLength=0.0125;
% ylabel(hcb,'\DeltaAbs (mOD)','FontWeight','bold','FontSize',12);

% Label other axes
xlabel(['Wavenumbers' ' (cm^{–1})'],'FontWeight','bold')
ylabel('Delays (ns)','FontWeight','bold')
hline = refline(0,0); hline.Color = [0.5 0.5 0.5];
switch TitleOnOff
    case 'On'
        title({datafilename;[rawcorr,' DATA';'']},'Interpreter','none')  
end
switch TitleOnOff
    case 'On'
        title({datafilename;[rawcorr,' DATA';'']},'Interpreter','none')  
end
