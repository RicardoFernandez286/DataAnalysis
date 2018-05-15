
%% READ inputs
% datedir     = 'D:\Data\20180314';
% samplename  = 'H2O_FeCN_dil_225817';
% blankname   = 'H2O_230247';
clear all

datedir     = 'D:\Data\20180511';
samplename  = 'R06_mTiO2_Re1213_MeOH_IRinsitu_150752';
blankname   = 'blank_FTIR_pyro_143537';

filt_order  = 8;


if exist([datedir filesep 'CalibratedProbe.csv'],'file') ~= 0
    wavenumbers = csvread([datedir filesep 'CalibratedProbe.csv']);
else
    wavenumbers = csvread([datedir filesep samplename filesep samplename '_wavenumbers.csv']);
end

sample_spec(:,1) = csvread([datedir filesep samplename filesep samplename '_probe_LinSpec.csv']);
sample_spec(:,2) = csvread([datedir filesep samplename filesep samplename '_reference_LinSpec.csv']);

corr_probe  = sample_spec(:,1) - medfilt1(sample_spec(:,1),filt_order,[],'truncate');
corr_ref    = sample_spec(:,2) - medfilt1(sample_spec(:,2),filt_order,[],'truncate');

p1 = plot(wavenumbers,corr_probe-min(corr_probe),'-','Color',[0 1 0],'DisplayName','Probe');
hold on
p2 = plot(wavenumbers,corr_ref-min(corr_ref),'-','Color',[1 0 0],'DisplayName','Reference');
hold off

rl = refline(0,0);
rl.Color = [0.5 0.5 0.5];

legend([p1,p2])
legend('boxoff')

xlabel('Wavenumbers (cm^{-1})','FontWeight','bold');
ylabel('Absorbance (mOD)','FontWeight','bold');