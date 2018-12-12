%% Read the data and calculate basic parameters
% Settings
Ncontours = 50;

% Read the data
folder      = '\\idnetapp-chem.uzh.ch\g_chem_hamm$\Group\Data from instruments\Lab 2 - TRUVIS\2018\Probe tests';
dataset     = 'SapphireWL_1000Shots_203548';

FNstart     = [folder filesep dataset filesep dataset];

wavelengths = csvread([FNstart '_wavelengths.csv']);
probe_shots = csvread([FNstart '_sp0_ds0_intensity_0.csv'])';
ref_shots   = csvread([FNstart '_sp0_ds1_intensity_0.csv'])';

% Calculate <I>
mean_probe  = mean(probe_shots,1);
mean_ref    = mean(ref_shots,1);

% Calculate δI
dIprobe     = probe_shots - mean(probe_shots,1);
dIref       = ref_shots - mean(ref_shots,1);

%% Calculate the correlation matrices (for δI)
dPP_Corr = corr(dIprobe,dIprobe);
dPR_Corr = corr(dIprobe,dIref);

figure
contourf(dPP_Corr,Ncontours,'EdgeColor','flat');
colorbar
caxis([-1 1])
pbaspect([1 1 1])
colormap(darkb2r(-1,1,Ncontours,2))
dline = refline(1,0);
dline.Color = 'w';
xlabel('Probe array pixel');
ylabel('Probe array pixel');

figure
contourf(dPR_Corr,Ncontours,'EdgeColor','flat');
colorbar
caxis([-1 1])
pbaspect([1 1 1])
colormap(othercolor('BuDRd_18',)
dline = refline(1,0);
dline.Color = 'w';
xlabel('Probe array pixel');
ylabel('Reference array pixel');


%% Init figure
fh = figure;
ax = axes(fh);
fh.Units    = 'normalized';
fh.Position = [0.25 0.4 0.6 0.5];

%% Get one pixel
Select_WL   = 468; 
pix         = findClosestId2Val(wavelengths,Select_WL);

%% Show trajectory for that pixel
reprate             = 2500; % Hz
Nshots              = size(probe_shots,1);
time                = (0:Nshots-1)./reprate; % in seconds

ax=subplot(2,4,1,ax);
plot(ax,time,dIprobe(:,pix),'g')
ax=subplot(2,4,5);
plot(ax,time,dIref(:,pix),'r')

%% Show histogram for that pixel
Nbins               = 50;

ax=subplot(2,4,2);
histogram(ax,dIprobe(:,pix),Nbins,'FaceColor','g')
xl = xlim;
xlim([-max(abs(xl)),max(abs(xl))])

ax=subplot(2,4,6);
histogram(ax,dIref(:,pix),Nbins,'FaceColor','r')
xl = xlim;
xlim([-max(abs(xl)),max(abs(xl))])

%% Show autocorrelation for that pixel
Nlags = Nshots-1;
[acf_dIprobe,lags]    = autocorr(dIprobe(:,pix),Nlags);
acf_dIref             = autocorr(dIref(:,pix),Nlags);
lags=lags./reprate;
plot(lags,acf_dIprobe)
NoiseSpec_Probe = abs(fft(acf_dIprobe));
loglog(NoiseSpec_Probe(1:round(Nshots/2)))

%% Show correlation coefficient between probe and reference arrays
Nlags           = Nshots/2-1;
DeltaIprobe     = downsample(dIprobe,2,0)-downsample(dIprobe,2,1);
DeltaIref       = downsample(dIref,2,0)-downsample(dIref,2,1);
ax=gca;
histogram(ax,DeltaIprobe(:,pix),Nbins,'FaceColor','g')
xl = xlim;
xlim([-max(abs(xl)),max(abs(xl))])

% [acf_DeltaIprobe,lags] = autocorr(DeltaIprobe(:,pix),Nlags);
% plot(lags,acf_DeltaIprobe)
% loglog(abs(fft(acf_DeltaIprobe)))

DeltaPP_Corr = corr(DeltaIprobe,DeltaIprobe);
DeltaPR_Corr = corr(DeltaIprobe,DeltaIref);
% 
% contourf(PR_Corr); colorbar
% dline = refline(1,0);
% dline.Color = 'w';

