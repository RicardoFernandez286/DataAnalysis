clear all
tic

c_0 = 2.99792458e-5;    % Speed of light in cm/fs
%%
rootdir = 'D:\GoogleDrive\RESULTS\From James\2DIR - 16_[Mn(bpy-ester)(CO)3Br]_Sample2_DCM_TD4FrPhCCl_14fs_20s';
datafilename = 'IRPS_IR_20210318_182451';

AllData = readmatrix([rootdir filesep datafilename '.2D'],'FileType','delimitedtext');  % read all data
t2delays = readmatrix([rootdir filesep datafilename '.DT'],'FileType','delimitedtext'); % in ps

%%
Ndelays = length(t2delays);
Npixels = size(AllData,2);
Nbins   = size(AllData,1)/Ndelays;

dt1 = 14; % fs
%%
TD_data = zeros(Nbins,Npixels,Ndelays);
ApoFunc = cos(pi*(1:Nbins)/(2*Nbins))';

for i=1:Ndelays
    TD_data(:,:,i) = ApoFunc.*AllData((Nbins*(i-1)+1):(Nbins*i),:);
end

%%
proberng = [1:80 100:Npixels];
delay = 25;
w0 = 1670;
zeropad_factor = 2;

% N_FTpoints = zeropad_factor*Nbins;

N_FTpoints = 1024;

ProbeAxis = 1:Npixels; ProbeAxis = ProbeAxis(proberng);

Resolution = 1/(N_FTpoints*dt1*c_0);
PumpAxis   = ((1:1:N_FTpoints)-1)'.*Resolution+w0;
% w1 = 1:Nbins * 

proc2D = real(fft(TD_data(:,proberng,delay),N_FTpoints));

contourf(PumpAxis,ProbeAxis,proc2D',40,'EdgeColor','flat');
colorbar

%%
toc