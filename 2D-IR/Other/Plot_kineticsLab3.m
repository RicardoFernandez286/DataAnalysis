%% Hardcoded settings and options
clear all
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

DoPlot      = 2; % 0 = plot nothing, 1 = plot everything, 2= plot what you select
whichplot   = 27;
fit_type    = 'No exp.'; % 'Single exp.' or 'Double exp.'
m           = 6;
n           = 3;

% Physical constants
HeNe        = 2.11079;          % HeNe period (fs)
c_0         = 2.99792458e-5;    % Speed of light in cm/fs


%% Get the list of files to be read
datafolder  = uigetdir();

% Get the list of FOLDERS in the folder
filelist    = dir(datafolder);
filenames   = {filelist.name}';
filenames   = filenames([filelist.isdir]);
filenames   = filenames(3:end);
filenames   = filenames(~contains(filenames,'Scattering'));

% Read the formatted strings of each filename
L           = length(filenames);
time        = zeros(1,L);
for i=1:L
    parts       = strsplit(filenames{i},'_');
    time_str    = strsplit(parts{end-2},'fs');
    time(i)     = str2double(time_str{1})./1000;
end

% Load files [scattering]
NR_scatt_fn  = dir([datafolder filesep 'Scattering' filesep '*NR_TD*']);
RE_scatt_fn  = dir([datafolder filesep 'Scattering' filesep '*RE_TD*']);
NR_scatt     = csvread([NR_scatt_fn.folder filesep NR_scatt_fn.name])';
RE_scatt     = csvread([RE_scatt_fn.folder filesep RE_scatt_fn.name])';

%% Read the pinhole data and do the W3 calibration
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

% Fit w3 axis
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

%% Read and process the 2D data
for i=1:L
    % Read the data and store them in a matrix: W1xW3xt2
    datafile_RE         = dir([datafolder filesep filenames{i} filesep '*RE_TD*']);
    datafile_NR         = dir([datafolder filesep filenames{i} filesep '*NR_TD*']);
    tempdata_RE(:,:,i)  = (csvread([datafile_RE.folder filesep datafile_RE.name]))';
    tempdata_NR(:,:,i)  = (csvread([datafile_NR.folder filesep datafile_NR.name]))';
    data_RE(:,:,i)      = tempdata_RE(:,1:end-1,i);
    data_NR(:,:,i)      = tempdata_RE(:,1:end-1,i);
    % Process the data
    RE_FT(:,:,i)        = fft(data_RE(:,:,i),0.5*Nbins*zeropad_factor,1);
    NR_FT(:,:,i)        = fft(data_NR(:,:,i),0.5*Nbins*zeropad_factor,1);
    Magnitude2D(:,:,i)  = abs(RE_FT(:,:,i) + NR_FT(:,:,i));
    % Phased2D(:,:,i)     = real((RE_FT(:,:,i) + NR_FT(:,:,i)).*exp(-1i*Phase_angle));  ?????
end

% Calculate the pump axis
PumpAxis = ((1:1:Nbins*0.5*zeropad_factor)-1)./(Nbins*0.5*zeropad_factor*HeNe*c_0);

% Define the cut region
cut_limits  = findClosestId2Val(PumpAxis,[ProbeAxis(1) ProbeAxis(end)]);

%% Plot
if DoPlot == 1
    for i=1:L
        % Define the axis
        ax = subplot(m,n,i);
        % Define limits and contour levels
        min_cut         = min(min(Magnitude2D(cut_limits(1):cut_limits(2),:,i)))*plot_percentcolrg/100;
        max_cut         = max(max(Magnitude2D(cut_limits(1):cut_limits(2),:,i)))*plot_percentcolrg/100;

        if plot_symmcolrange
            min_cut     = -max(abs([min_cut max_cut]));
            max_cut     = max(abs([min_cut max_cut]));
        end

        Contours        = linspace(min_cut,max_cut,Ncontours);

        % Make the contour plot
        contourf(ax,PumpAxis,ProbeAxis,Magnitude2D(:,:,i)',Contours,'LineColor','flat')

        % Show the contour lines every Xth level
        if plot_showlines
            % Define the contours to plot
            step = ((max_cut - min_cut)./Ncontours).*plot_skiplevels;
            neg_zero = (plot_Nwhites/Ncontours)*min_cut;
            pos_zero = (plot_Nwhites/Ncontours)*max_cut;

            plot_contours = [min_cut:step:neg_zero pos_zero:step:max_cut];

            hold(ax,'on')
            contour(ax,PumpAxis,ProbeAxis,Magnitude2D(:,:,i)',plot_contours,'LineColor',0.2*[1 1 1],'LineStyle','-','LineWidth',0.1);
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

    %     % Add phase text
    %     text(ax,50,300,['t2 ' num2str(time(i)) ' ps'],'Units','pixels','FontSize',12);
    end
elseif DoPlot == 2
        i = whichplot;
        ax=gca;
        % Define limits and contour levels
        min_cut         = min(min(Magnitude2D(cut_limits(1):cut_limits(2),:,i)))*plot_percentcolrg/100;
        max_cut         = max(max(Magnitude2D(cut_limits(1):cut_limits(2),:,i)))*plot_percentcolrg/100;

        if plot_symmcolrange
            min_cut     = -max(abs([min_cut max_cut]));
            max_cut     = max(abs([min_cut max_cut]));
        end

        Contours        = linspace(min_cut,max_cut,Ncontours);

        % Make the contour plot
        contourf(ax,PumpAxis,ProbeAxis,Magnitude2D(:,:,i)',Contours,'LineColor','flat')

        % Show the contour lines every Xth level
        if plot_showlines
            % Define the contours to plot
            step = ((max_cut - min_cut)./Ncontours).*plot_skiplevels;
            neg_zero = (plot_Nwhites/Ncontours)*min_cut;
            pos_zero = (plot_Nwhites/Ncontours)*max_cut;

            plot_contours = [min_cut:step:neg_zero pos_zero:step:max_cut];

            hold(ax,'on')
            contour(ax,PumpAxis,ProbeAxis,Magnitude2D(:,:,i)',plot_contours,'LineColor',0.2*[1 1 1],'LineStyle','-','LineWidth',0.1);
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

    %     % Add phase text
    %     text(ax,50,300,['t2 ' num2str(time(i)) ' ps'],'Units','pixels','FontSize',12);
end

%% Plot kinetics

SelTraces = inputdlg('Enter the coordinates of the points to plot along t2 (format: x1,y1;x2,y2...):',...
                 'Input desired coordinates', [1 60]);
if isempty(SelTraces)
    return
end
SelTraces       = str2num(SelTraces{:});
Q               = size(SelTraces,1);
pump_search     = SelTraces(:,1); 
probe_search    = SelTraces(:,2);

for i=1:Q
    probe_index(i)      = findClosestId2Val(ProbeAxis,probe_search(i));
    pump_index(i)       = findClosestId2Val(PumpAxis,pump_search(i));
end

% Get the data
for i=1:Q
    kindata(:,i) = Magnitude2D(pump_index(i),probe_index(i),:);
    caption{i} = ['(' num2str(round(PumpAxis(pump_index(i)))) ', ' num2str(round(ProbeAxis(probe_index(i)))) ') cm^{-1}'];
end

% Create a new figure with consistent format
    label = '2D signal (a.u.)';
    fh = figure();
    fh.Position(3)  = 800;
    fh.Position(4)  = 425;
    fh.Color        = [1 1 1];
    % Define the axes
    axes2 = axes('Parent',fh);
    axes(axes2);
    cmap=colormap(othercolor('Mrainbow',Q));
    % Plot the data
    for n=1:Q
       plot(axes2,time,kindata(:,n),'o','LineWidth',2,'MarkerSize',2,'color',cmap(n,:));
       hold on
    end
    %%% Nice formatting
    set(gca,'FontSize',14)
    xlabel('t_{2} delay (ps)','FontSize',14,'FontWeight','bold');
    ylabel(label,'FontSize',14,'FontWeight','bold')
    % title([handles.datafilename;'WAITING TIME KINETICS';''],'Interpreter','none','FontSize',10)
    axis tight

    % Show only positive t2 times
    xlim([0 max(time)]);
    
    % Show a bit below zero
    ylim([-0.1.*max(kindata(:)) max(kindata(:))])

    % Create legend
    legend(gca,caption,'FontSize',12)
    legend('boxoff')
    legend('Location','northeast')

    % Set linear or log X scale
    set(axes2,'xscale','lin');

    % Add zero line
    hline = refline(axes2,[0 0]); hline.Color = [0.5 0.5 0.5];
    set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

    axes2.Units     = 'pixels';
    axes2.Position  = [75 75 675 320];
    axes2.Units     = 'normalized';
    
%% Fit the data
for i=1:Q
    switch fit_type
        case 'Double exp.'
            fit_func   = @(P,t) P(1)*exp(-t/P(2))+P(3)*exp(-t/P(4));
            LB = [0 0 0 0];
            UB = [Inf Inf Inf Inf];
            start_param = [600 10 600 10];
        case 'Single exp.'
            fit_func   = @(P,t) P(1)*exp(-t/P(2));
            LB = [0 0];
            UB = [Inf Inf];
            start_param = [600 2];
    end
    options = optimoptions('lsqcurvefit',...
                'MaxFunctionEvaluations',1000,...
                'MaxIterations',1000,...
                'Algorithm','trust-region-reflective',...  %'levenberg-marquardt' 'trust-region-reflective'
                'OptimalityTolerance',1e-15,...
                'FunctionTolerance',1e-15,...
                'StepTolerance',1e-15);
            
    [fitparam(i,:),~,residuals,~,~,~,jacobian] = lsqcurvefit(fit_func,start_param,time(time>=0),kindata(time>=0,i)',LB,UB,options);
    fitted_CI           = nlparci(fitparam(i,:),residuals,'jacobian',jacobian);
    fitted_err(i,:)     = ((fitted_CI(:,2) - fitted_CI(:,1))/2)';
end

for n=1:Q
    plot_time = linspace(0,max(time),1000);
    plot(axes2,plot_time,fit_func(fitparam(n,:),plot_time),'-','LineWidth',1,'color',cmap(n,:));
    hold on
    disp('==================================')
    disp(['Fit results for Curve ' num2str(n)])
    disp(['Tau1 = ' num2str(fitparam(n,2)) ' ± ' num2str(fitted_err(n,2)) ' ps'])
    disp(['A1   = ' num2str(fitparam(n,1)) ' ± ' num2str(fitted_err(n,1))])
    switch fit_type
        case 'Double exp.'
        	disp(['Tau2 = ' num2str(fitparam(n,4)) ' ± ' num2str(fitted_err(n,4)) ' ps'])
            disp(['A2   = ' num2str(fitparam(n,3)) ' ± ' num2str(fitted_err(n,3))])
    end
    disp('==================================')
end
hold off

csvwrite('Kinetic_data.csv',[time' kindata])