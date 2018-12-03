function Test_lab3_phasing(ax,Phase_angle)
% Let's get estarted.
% Phase_angle = pi/1;

%% Initialization
% Hardcoded settings
probe_fitorder  = 2;    % Order of the polynomial to fit to the scattering calibration
zeropad_factor  = 1;    % Zeropad the data by a factor X
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
datafolder  = '\\idnetapp-chem.uzh.ch\g_chem_hamm$\Group\Data from instruments\Lab 3 - 2DIR\20181022\Parallel1752_ce_LiNO3_KSCN2iso_2mM_1750fs_150_2025';

% Load files [Sample]
NR_data    = csvread([datafolder filesep 'Parallel1823_ce_LiNO3_KSCN2iso_2mM_7500fs_150_2025_NR_TD_20_spectra_.csv'])';
RE_data    = csvread([datafolder filesep 'Parallel1752_ce_LiNO3_KSCN2iso_2mM_1750fs_150_2025_RE_TD_27_spectra_.csv'])';

% Load files [scattering]
NR_scatt     = csvread([datafolder filesep 'Parallel1700_pinhole_300fs_150_2025_NR_TD_7_spectra_.csv'])';
RE_scatt     = csvread([datafolder filesep 'Parallel1700_pinhole_300fs_150_2025_RE_TD_8_spectra_.csv'])';

%% Preprocess
% Define t1 time axis
Nbins       = size(RE_scatt,1) + size(NR_scatt,1);
t1_delays   = linspace(-Nbins/2,Nbins/2,Nbins).*HeNe; % in fs

% Remove "pixel 33"
RE_scatt    = RE_scatt(:,1:end-1);
NR_scatt    = NR_scatt(:,1:end-1);

% Average values at time zero
    % I have no idea

% Join the NR and RE arrays
NRRE_scatt  = [NR_scatt;RE_scatt];

% Calculate w1 axis
PumpAxis = ((1:1:Nbins)-1)./(Nbins*HeNe*c_0);

% Calculate the power spectrum
Power_spec = abs(fft(NRRE_scatt));

%% Fit w3 axis
[~,binspecmax]          = max(sum(Power_spec,2));
P                       = floor(Nbins/50);
fitrange                = (binspecmax-P):(binspecmax+P);
Npixels                 = size(NRRE_scatt,2);
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

% Remove "pixel 33"
NR_data     = NR_data(:,1:end-1);
RE_data     = RE_data(:,1:end-1);

% Join the two arrays
NRRE_data   = [NR_data;RE_data];

% Fourier transform
NR_FT       = fft(NR_data,0.5*Nbins*zeropad_factor,1);
RE_FT       = fft(RE_data,0.5*Nbins*zeropad_factor,1);

% Calculate the pump axis
PumpAxis = ((1:1:Nbins*0.5*zeropad_factor)-1)./(Nbins*0.5*zeropad_factor*HeNe*c_0);

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