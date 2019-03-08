function Test_lab3_phasing(ax,Phase_angle)
% Let's get estarted.
% Phase_angle = pi/1;

%% Initialization
% Hardcoded settings
probe_fitorder  = 2;    % Order of the polynomial to fit to the scattering calibration
zeropad_factor  = 1;    % Zeropad the data by a factor X
cos_exp         = 1;
Ncontours       = 50;   % No. of contours to plot
plot_Nwhites    = 2;    % No. of white levels
plot_showlines  = 0;    % Show black lines or not (1 or 0)
plot_skiplevels = 2;    % No. of levels to skip for plotting the black lines   
plot_symmcolrange=1;    % Plot with a symmetric colour range (1 or 0) - default = 1
plot_percentcolrg=100;

% Physical constants
HeNe        = 2.11079;          % HeNe period (fs)
c_0         = 2.99792458e-5;    % Speed of light in cm/fs

%% Load files
% Define the data folder
rootdir     = '\\idnetapp-chem.uzh.ch\g_chem_hamm$\Group\Data from instruments\Lab 3 - 2DIR\20181216\TwoDSpectra';
dataname    = 'Parallel2244_solution_T19_75000fs_150_2025_TD';

% Do the file reading
filelist    = dir(rootdir);
filenames   = {filelist.name}';
filenames   = filenames(~[filelist.isdir]);
pinholeFile = [rootdir filesep filenames{contains(filenames,'Pinhole','IgnoreCase',true)}];
filename    = [rootdir filesep dataname '.csv'];

% Open the files
ds1         = importdata(filename,',',1);
ds2         = importdata(filename,',',36);
ds3         = importdata(filename,',',71);
ds4         = importdata(filename,',',106);

ds1_pinhole = importdata(pinholeFile,',',1);
ds2_pinhole = importdata(pinholeFile,',',36);
ds3_pinhole = importdata(pinholeFile,',',71);
ds4_pinhole = importdata(pinholeFile,',',106);

% Join the four states in a 3D array
data(:,:,1) = ds1.data(1:end-1,:)';
data(:,:,2) = ds2.data(1:end-1,:)';
data(:,:,3) = ds3.data(1:end-1,:)';
data(:,:,4) = ds4.data(1:end-1,:)';

pinhole(:,:,1) = ds1_pinhole.data(1:end-1,:)';
pinhole(:,:,2) = ds2_pinhole.data(1:end-1,:)';
pinhole(:,:,3) = ds3_pinhole.data(1:end-1,:)';
pinhole(:,:,4) = ds4_pinhole.data(1:end-1,:)';

% Initialize some variables
Nbins       = 0.5*(size(data,1)+1);
Npixels     = size(data,2);
t1delays    = linspace(0,Nbins-1,Nbins)*HeNe;

% Split the data into rephasing and non-rephasing
data_NR     = data(1:Nbins,:,:);
data_RE     = data(Nbins:end,:,:);
pinhole_NR  = pinhole(1:Nbins,:,:);
pinhole_RE  = pinhole(Nbins:end,:,:);

% Add the datastates to obtain the correct signal
% CORRECT ADDITION: 1 - 2 + 3 - 4
% signal_xx = Nbins x Npixels
signal_NR   = data_NR(:,:,1) - data_NR(:,:,2) + data_NR(:,:,3) - data_NR(:,:,4);
signal_RE   = data_RE(:,:,1) - data_RE(:,:,2) + data_RE(:,:,3) - data_RE(:,:,4);

pinhole_NR  = pinhole_NR(:,:,1) - pinhole_NR(:,:,2) + pinhole_NR(:,:,3) - pinhole_NR(:,:,4);
pinhole_RE  = pinhole_RE(:,:,1) - pinhole_RE(:,:,2) + pinhole_RE(:,:,3) - pinhole_RE(:,:,4);
% Calculate w1 axis
PumpAxis = ((1:1:Nbins)-1)./(Nbins*HeNe*c_0);

% Calculate the power spectrum
Power_spec = abs(fft(pinhole_NR))+abs(fft(pinhole_RE));

%% Fit w3 axis
[~,binspecmax]          = max(sum(Power_spec,2));
P                       = floor(Nbins/50);
fitrange                = (binspecmax-P):(binspecmax+P);
Npixels                 = size(signal_NR,2);
pixels                  = transpose(1:1:Npixels);
% Get the maxima for probe calibration    
[~,maxindex]            = max(Power_spec,[],1);
scattering_maxima       = PumpAxis(maxindex)';
switch probe_fitorder
    case 1
        model       = 'poly1';
    case 2
        model       = 'poly2';
end
warning('off','curvefit:fit:iterationLimitReached');
mdl                     = fit(pixels,scattering_maxima,model,'Robust','Bisquare');
ProbeAxis               = mdl(pixels);


%% Process the SAMPLE

% Calculate the pump axis
PumpAxis = ((1:1:Nbins*zeropad_factor)-1)./(Nbins*zeropad_factor*HeNe*c_0);
q=1:Nbins;
cosine             = ((cos(pi*(q)/(2*Nbins))).^cos_exp)';
box                = (heaviside(q))';
apodize_function   = cosine.*box;
apo_NR             = signal_NR.*apodize_function;
apo_RE             = signal_RE.*apodize_function;
% FFT and phase
NR_FT              = fft(apo_NR,Nbins*zeropad_factor);
RE_FT              = fft(apo_RE,Nbins*zeropad_factor);

% Phase the spectrum
Absorptive  = abs((NR_FT + RE_FT).*exp(-1i*Phase_angle));

% Define the cut region
cut_limits  = findClosestId2Val(PumpAxis,[ProbeAxis(1) ProbeAxis(end)]);

%% Plot
% Define limits and contour levels
min_cut    = min(min(Absorptive(cut_limits(1):cut_limits(2),:)))*plot_percentcolrg/100;
max_cut    = max(max(Absorptive(cut_limits(1):cut_limits(2),:)))*plot_percentcolrg/100;

if plot_symmcolrange
    min_cut = -max(abs([min_cut max_cut]));
    max_cut = max(abs([min_cut max_cut]));
end

Contours    = linspace(min_cut,max_cut,Ncontours);
    
% Make the contour plot
contourf(ax,PumpAxis,ProbeAxis,Absorptive',Contours,'LineColor','flat')

% Show the contour lines every Xth level
if plot_showlines
    % Define the contours to plot
    step = ((max_cut - min_cut)./Ncontours).*plot_skiplevels;
	neg_zero = (plot_Nwhites/Ncontours)*min_cut;
    pos_zero = (plot_Nwhites/Ncontours)*max_cut;

    plot_contours = [min_cut:step:neg_zero pos_zero:step:max_cut];

    hold(ax,'on')
    contour(ax,PumpAxis,ProbeAxis,Absorptive',plot_contours,'LineColor',0.2*[1 1 1],'LineStyle','-','LineWidth',0.1);
    hold(ax,'off')
end   

% Make uniform, consistent format
    ax.FontSize = 10;
    ax.LineWidth = 1;
    ax.TickLength = [0.015 0.035];
    ax.Visible='On';
    
% Colormap and colorbar
cmap    = darkb2r(min_cut,max_cut,Ncontours,plot_Nwhites);
colormap(cmap)
colorbar

% Axis titles
ylabel('\omega_3 (cm^{-1})','FontWeight','bold');
xlabel('\omega_1 (cm^{-1})','FontWeight','bold');

% Diagonal line
hline = refline(1,1);
hline.Color = [0 0 0];

% Axis limits
xlim([ProbeAxis(1) ProbeAxis(end)])
ylim([ProbeAxis(1) ProbeAxis(end)])

% Add phase text
% text(ax,50,300,['Phase: ' num2str(Phase_angle) ' rad'],'Units','pixels','FontSize',12);