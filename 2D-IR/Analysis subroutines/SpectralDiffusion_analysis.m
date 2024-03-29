function  dataStruct = SpectralDiffusion_analysis(app,dataStruct,Axes,k,varargin)

% Description: This function dataStruct the spectral diffusion analysis of a
% series of 2D-IR data as a function of the waiting time
% 
% Usage: dataStruct = SpectralDiffusion_analysis(dataStruct)
% Inputs:
%     datatype ('Raw' or 'Signal', as in the MESS program)
%     datafilename
%     rootdir
%     
% Outputs:
%     cmprobe           (Double)
%     bins              (Double)
%     t2delays          (Double)
%     Ndelays           (Double)
%     Nspectra          (Double)
%     Ndatastates       (Double)
%     Nbins             (Double)
%     Nslowmod          (Double)
%     count             (Cell array)
%     probe             (Cell array)
%     reference         (Cell array)
%     interferogram     (Cell array)
%     signal			(Cell array)
%
% Ricardo Fern�ndez-Ter�n / 15.04.2019 / v2.0a

%% READ from dataStruct
% Get the data
ProbeAxis           = dataStruct.ProbeAxis{k};
PumpAxis            = dataStruct.PumpAxis;
t2delays            = dataStruct.t2delays;
Ndelays             = length(t2delays);

if ~isempty(varargin)
    mode        = varargin{1};
    settings    = varargin{2};
    saveTraces  = 'off';
    t2_fitrange = settings.t2_fitrange;
    
    t2_lim      = sort(findClosestId2Val(t2delays,t2_fitrange));
    t2_idx      = t2_lim(1):1:t2_lim(2);
    Ndelays_anl = length(t2_idx);
else
    mode = 'Single';
    % Get information from the GUI
    plot_pumpdirection  = app.I2D_PumpAxisOrientation.Value;
    saveTraces          = app.I2D_SavetracesSwitch.Value;
    interactivemode     = app.I2D_InteractivemodeSwitch.Value;
    
    fitrange = inputdlg('Input t2 range to fit (format: start,end)','Select t2 range',[1 50]);
    t2_fitrange = str2num(fitrange{:});
    
    t2_lim      = sort(findClosestId2Val(t2delays,t2_fitrange));
    t2_idx      = t2_lim(1):1:t2_lim(2);
    Ndelays_anl = length(t2_idx);
end


PROC_2D_DATA        = dataStruct.PROC_2D_DATA;
Npixels             = length(ProbeAxis);

% Hardcoded settings
Interp_method       = 'InterpFT'; % 2D FFT, InterpFT, Mesh, None. Default = InterpFT
% Fit_type            = 'Quadratic'; % 'Quadratic' 'None' - Sets the type of fit to calculate the minima along the slices

Interp_order        = 6; % Multiplies Npixels by a factor to interpolate in the probe dimension
Nfitpoints          = 4;
N_pointsParabola    = round(Interp_order*Nfitpoints); % No. of X points to fit a parabola on each side of the peak: [-x (peak) +x]
N_points_IvCLS      = 8;

intensity_threshold = 60/100; % Intensity threshold for the linear fit to get the CLS or IvCLS
textcolor           = 'none'; % 'none' or RGB color
fit_method          = 'LSQ'; % 'GA' or 'LSQ' for Genetic algorithm or Least-squares

diagnosticPlot      = 0; % Plots the cuts and fits for diagnostic purposes

% Suppress warnings
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');

% Initialize variables (ensures only the current analysis is saved)
CLS_Xdata     = {};
CLS_Ydata     = {};
IvCLS_Xdata   = {};
IvCLS_Ydata   = {};
NLS_Xdata     = {};
NLS_Ydata     = {};

%% Ask the user which method to use (if not in 'Batch' mode)
specdif_options     = {'CLS','IvCLS','NLS','CLS+IvCLS','Gaussian Fit'};

if strcmp(mode,'Batch')
    doAnalysis       = 1;
    specdif_typeindx = 1;
    
    diagnosticPlot   = 0;
    fit_method       = 'None';
    
    intensity_threshold = settings.intensity_threshold;
    Interp_method       = settings.Interp_method;
    interactivemode     = 'Batch';
else
    [specdif_typeindx,doAnalysis] = listdlg('ListString',specdif_options,'OKstring','Continue','SelectionMode','single','ListSize',[150,100],'PromptString','Select analysis type:');
end

if doAnalysis == 0
    return
end

%% Interpolate the data along the probe axis (interpolation along w1 can be performed by apodisation)
Interp_Data = cell(Ndelays,1);
switch Interp_method
    case '2D FFT'
        for m=1:Ndelays
            Interp_Data{m,k}    = real(fft(ifft(PROC_2D_DATA{m,k},2),Npixels.*Interp_order,2));
        end
    case 'InterpFT'
        for m=1:Ndelays
            Interp_Data{m,k}    = interpft(PROC_2D_DATA{m,k},Npixels.*Interp_order,2);
        end
    case 'Mesh'
        % Do something
    case 'None'
        for m=1:Ndelays
            Interp_Data{m,k} = PROC_2D_DATA{m,k};
        end
end

%% Get the spectral diffusion dynamics
% Interpolate the probe axis
pixels              = 1:1:Npixels;
Newpixels           = linspace(1,Npixels,Npixels.*Interp_order);
Interp_ProbeAxis    = interp1(pixels,ProbeAxis,Newpixels);

switch interactivemode
    case 'On'
        % Select points interactively, to finish hit RETURN
        [RegionPositions,~] = SelectRegions(Axes);
        % Separate the values into pump and probe, according to how the data is plotted
        % i.e. convert from [x_min x_max y_min y_max] to [pump_min pump_max probe_min probe_max]
        switch plot_pumpdirection
            case 'Horizontal'
                pump_ranges     = RegionPositions(:,1:2); 
                probe_ranges    = RegionPositions(:,3:4);
            case 'Vertical'
                pump_ranges     = RegionPositions(:,3:4); 
                probe_ranges    = RegionPositions(:,1:2);
        end
    case 'Off'
        RegionPositions_cell    = inputdlg({'Enter the pump regions (Format: min1 max1; min2 max2; ...): ','Enter the probe regions (Format: min1 max1; min2 max2; ...):'},'Define regions to integrate', [1 80]);
        pump_ranges             = str2num(RegionPositions_cell{1}); %#ok<*ST2NM>
        probe_ranges            = str2num(RegionPositions_cell{2});
    case 'Batch'
        pump_ranges             = settings.pump_ranges;
        probe_ranges            = settings.probe_ranges;
end

% Check that we have the correct sizes
L = size(pump_ranges,1);
M = size(probe_ranges,1);
if L ~= M
    warndlg('Number of pump and probe ranges is different!')
    return
end

SpecDif_ind_all = zeros(Ndelays,L);

for q=1:L
    % Parse the pump ranges into indexes
    pump_idxrange   = findClosestId2Val(PumpAxis{1,1},pump_ranges(q,:));
    probe_idxrange  = findClosestId2Val(Interp_ProbeAxis,probe_ranges(q,:));

    % Prepare figure and axes for the diagnostic plots
    if diagnosticPlot == 1
        diagnosticFig   = figure;
        diagnosticAx    = axes('parent',diagnosticFig);
        cla(diagnosticAx,'reset');
        hold(diagnosticAx,'on')
    end

    % Get the data according to the selected method
    switch specdif_options{specdif_typeindx}
        case 'CLS'
            % Minima in w3 for each value of w1
            probe_indexes       = probe_idxrange(1):1:probe_idxrange(end);
            pump_indexes        = pump_idxrange(1):1:pump_idxrange(end);
            Npeaks              = length(pump_indexes);
            min_values          = NaN*zeros(Npeaks,Ndelays);
            SpecDif_ind         = NaN*zeros(length(t2delays),1);
            pump_cut            = cell(Ndelays,1);
            CLS_Xdata           = cell(Ndelays,1);
            CLS_Ydata           = cell(Ndelays,1);
            for m=t2_idx
                % Get the data around the minimum
                [peak_val,peak_pos] = min(Interp_Data{m,k}(pump_indexes,probe_indexes),[],2);
                peak_pos            = peak_pos + probe_idxrange(1);
                % Fit a quadratic polynomial within N points around the minimum 
                %   and find the analytical vertex of the parabola
                if diagnosticPlot == 1
                    cla(diagnosticAx,'reset'); hold(diagnosticAx,'on'); 
                end
                for i=1:Npeaks
                    start_probe     = max(peak_pos(i)-N_pointsParabola,1);
                    end_probe       = min(peak_pos(i)+N_pointsParabola,length(Interp_ProbeAxis));
                    probe_segment   = Interp_ProbeAxis(start_probe:end_probe);
                    data_segment    = Interp_Data{m,k}(pump_indexes(i),start_probe:end_probe);
                    coeff           = polyfit(probe_segment,data_segment,6);
                    fitData = polyval(coeff,probe_segment);
                    [~,min_idx]     = min(fitData);
                    min_values(i,m) = probe_segment(min_idx);
                    % min_values(i,m) = -coeff(2)/(2*coeff(1));
                    if diagnosticPlot == 1
                        cmap= jet(Npeaks);
                        plot(diagnosticAx,probe_segment,data_segment./max(abs(data_segment(:))),'o','MarkerSize',3,'Color',cmap(i,:));
                        plot(diagnosticAx,probe_segment,polyval(coeff,probe_segment)./max(abs(data_segment(:))),'-','LineWidth',2,'Color',cmap(i,:));
                        pause(0.05);
                    end
                end
                % Fit a line around the minimum according to the intensity cutoff
                peak_val            = peak_val/max(abs(peak_val));
                above_threshold     = abs(peak_val) >= max(abs(peak_val)).*intensity_threshold;
                pump_cut{m}         = pump_indexes(above_threshold);
                mdl                 = fit(PumpAxis{1,1}(pump_cut{m}),min_values(above_threshold,m),'poly1','Robust','Bisquare');
                CLS_coeff           = coeffvalues(mdl);
                SpecDif_ind(m)      = CLS_coeff(1);
                % Save the spectral diffusion data
                CLS_Xdata{m,1}        = PumpAxis{1,1}(pump_cut{m});
                CLS_Ydata{m,1}        = mdl(CLS_Xdata{m,1});
%                 CLS_Ydata{m,1}        = min_values(above_threshold,m);
            end
            SpecDif_name = 'CLS';
            % Save CLS to file
            if strcmp(saveTraces,'On')
                csvwrite([char(dataStruct.rootdir) filesep char(dataStruct.datafilename) '_CLS.csv'],[dataStruct.t2delays,SpecDif_ind])
            end

        case 'IvCLS'
            % Minima in w1 for each value of w3
            probe_indexes       = probe_idxrange(1):1:probe_idxrange(end);
            pump_indexes        = pump_idxrange(1):1:pump_idxrange(end);
            Npeaks              = length(probe_indexes);
            min_values          = zeros(Npeaks,Ndelays);
            SpecDif_ind         = zeros(Ndelays,1);
            probe_cut           = cell(Ndelays,1);
            IvCLS_Xdata         = cell(Ndelays,1);
            IvCLS_Ydata         = cell(Ndelays,1);
            for m=t2_idx
                % Get the data around the minimum
                [peak_val,peak_pos] = min(Interp_Data{m,k}(pump_indexes,probe_indexes),[],1);
                peak_pos            = peak_pos + pump_idxrange(1);
                % Fit a quadratic polynomial within 2 points around the minimum 
                %   and find the analytical vertex of the parabola
                if diagnosticPlot == 1
                    cla(diagnosticAx,'reset'); hold(diagnosticAx,'on'); 
                end
                for i=1:Npeaks
                    pump_segment    = PumpAxis{1,1}(peak_pos(i)-N_points_IvCLS:peak_pos(i)+N_points_IvCLS);
                    data_segment    = Interp_Data{m,k}(peak_pos(i)-N_points_IvCLS:peak_pos(i)+N_points_IvCLS,probe_indexes(i));
                    coeff           = polyfit(pump_segment,data_segment,2);
                    min_values(i,m) = -coeff(2)/(2*coeff(1));
                    if diagnosticPlot == 1
                        cmap= jet(Npeaks);
                        plot(diagnosticAx,pump_segment,data_segment./max(abs(data_segment(:))),'o','MarkerSize',3,'Color',cmap(i,:));
                        plot(diagnosticAx,pump_segment,polyval(coeff,pump_segment)./max(abs(data_segment(:))),'-','LineWidth',2,'Color',cmap(i,:));
                    end
                end
                % Fit a line around the minimum according to the intensity cutoff
                peak_val            = peak_val/max(abs(peak_val));
                above_threshold     = abs(peak_val) >= max(abs(peak_val)).*intensity_threshold;
                probe_cut{m}        = probe_indexes(above_threshold);
                mdl                 = fit(min_values(above_threshold,m),Interp_ProbeAxis(probe_cut{m})','poly1','Robust','Bisquare');
                IvCLS_coeff         = coeffvalues(mdl);
                SpecDif_ind(m)      = 1./IvCLS_coeff(1);
                % Save the spectral diffusion data
                IvCLS_Xdata{m}      = min_values(above_threshold,m);
                IvCLS_Ydata{m}      = Interp_ProbeAxis(probe_cut{m});
            end
            SpecDif_name    = 'IvCLS';
            % Save IvCLS to file
            if strcmp(saveTraces,'On')
                csvwrite([char(dataStruct.rootdir) filesep char(dataStruct.datafilename) '_IvCLS.csv'],[dataStruct.t2delays,SpecDif_ind])
            end

        case 'NLS'
            %%% Not yet implemented.

        case 'CLS+IvCLS'
            % IvCLS
            % Minima in w1 for each value of w3
            probe_indexes       = probe_idxrange(1):1:probe_idxrange(end);
            pump_indexes        = pump_idxrange(1):1:pump_idxrange(end);
            Npeaks              = length(probe_indexes);
            min_values          = zeros(Npeaks,Ndelays);
            IvCLS_value         = zeros(Ndelays,1);
            probe_cut           = cell(Ndelays,1);
            IvCLS_Xdata         = cell(Ndelays,1);
            IvCLS_Ydata         = cell(Ndelays,1);
            for m=t2_idx
                % Get the data around the minimum
                [peak_val,peak_pos] = min(Interp_Data{m,k}(pump_indexes,probe_indexes),[],1);
                peak_pos            = peak_pos + pump_idxrange(1);
                % Fit a quadratic polynomial within 2 points around the minimum 
                %   and find the analytical vertex of the parabola
                if diagnosticPlot == 1
                    cla(diagnosticAx,'reset'); hold(diagnosticAx,'on'); 
                end
                for i=1:Npeaks
                    pump_segment    = PumpAxis{1,1}(peak_pos(i)-N_points_IvCLS:peak_pos(i)+N_points_IvCLS);
                    data_segment    = Interp_Data{m,k}(peak_pos(i)-N_points_IvCLS:peak_pos(i)+N_points_IvCLS,probe_indexes(i));
                    coeff           = polyfit(pump_segment,data_segment,2);
                    min_values(i,m) = -coeff(2)/(2*coeff(1));
                    if diagnosticPlot == 1
                        cmap= jet(Npeaks);
                        plot(diagnosticAx,pump_segment,data_segment./max(abs(data_segment(:))),'o','MarkerSize',3,'Color',cmap(i,:));
                        plot(diagnosticAx,pump_segment,polyval(coeff,pump_segment)./max(abs(data_segment(:))),'-','LineWidth',2,'Color',cmap(i,:));
                    end
                end
                % Fit a line around the minimum according to the intensity cutoff
                peak_val            = peak_val/max(abs(peak_val));
                above_threshold     = abs(peak_val) >= max(abs(peak_val)).*intensity_threshold;
                probe_cut{m}        = probe_indexes(above_threshold);
                mdl                 = fit(min_values(above_threshold,m),Interp_ProbeAxis(probe_cut{m})','poly1','Robust','Bisquare');
                IvCLS_coeff         = coeffvalues(mdl);
                IvCLS_value(m)      = 1./IvCLS_coeff(1);
                % Save the spectral diffusion data
                IvCLS_Xdata{m}      = min_values(above_threshold,m);
                IvCLS_Ydata{m}      = Interp_ProbeAxis(probe_cut{m});
            end

            % CLS
            % Minima in w3 for each value of w1
            probe_indexes       = probe_idxrange(1):1:probe_idxrange(end);
            pump_indexes        = pump_idxrange(1):1:pump_idxrange(end);
            Npeaks              = length(pump_indexes);
            min_values          = zeros(Npeaks,Ndelays);
            CLS_value           = zeros(Ndelays,1);
            pump_cut            = cell(Ndelays,1);
            CLS_Xdata           = cell(Ndelays,1);
            CLS_Ydata           = cell(Ndelays,1);
            for m=t2_idx
                % Get the data around the minimum
                [peak_val,peak_pos] = min(Interp_Data{m,k}(pump_indexes,probe_indexes),[],2);
                peak_pos            = peak_pos + probe_idxrange(1);
                % Fit a quadratic polynomial within N points around the minimum 
                %   and find the analytical vertex of the parabola
                if diagnosticPlot == 1
                    cla(diagnosticAx,'reset'); hold(diagnosticAx,'on'); 
                end
                for i=1:Npeaks
                    probe_segment   = Interp_ProbeAxis(peak_pos(i)-N_pointsParabola:peak_pos(i)+N_pointsParabola);
                    data_segment    = Interp_Data{m,k}(pump_indexes(i),peak_pos(i)-N_pointsParabola:peak_pos(i)+N_pointsParabola);
                    coeff           = polyfit(probe_segment,data_segment,2);
                    min_values(i,m) = -coeff(2)/(2*coeff(1));
                    if diagnosticPlot == 1
                        cmap= jet(Npeaks);
                        plot(diagnosticAx,probe_segment,data_segment./max(abs(data_segment(:))),'o','MarkerSize',3,'Color',cmap(i,:));
                        plot(diagnosticAx,probe_segment,polyval(coeff,probe_segment)./max(abs(data_segment(:))),'-','LineWidth',2,'Color',cmap(i,:));
                    end
                end
                % Fit a line around the minimum according to the intensity cutoff
                peak_val            = peak_val/max(abs(peak_val));
                above_threshold     = abs(peak_val) >= max(abs(peak_val)).*intensity_threshold;
                pump_cut{m}         = pump_indexes(above_threshold);
                mdl                 = fit(PumpAxis{1,1}(pump_cut{m}),min_values(above_threshold,m),'poly1','Robust','Bisquare');
                CLS_coeff           = coeffvalues(mdl);
                CLS_value(m)        = CLS_coeff(1);
                % Save the spectral diffusion data
                CLS_Xdata{m}        = PumpAxis{1,1}(pump_cut{m});
                CLS_Ydata{m}        = min_values(above_threshold,m);
            end
            SpecDif_name = 'CLS';

            % Save CLS and IvCLS to file
            if strcmp(saveTraces,'On')
                csvwrite([char(dataStruct.rootdir) filesep char(dataStruct.datafilename) '_IvCLS.csv'],[dataStruct.t2delays,IvCLS_value])
                csvwrite([char(dataStruct.rootdir) filesep char(dataStruct.datafilename) '_CLS.csv'],[dataStruct.t2delays,CLS_value])
            end
    end

    SpecDif_ind_all(:,q) = SpecDif_ind;
end

if diagnosticPlot == 1
    hold(diagnosticAx,'off')
end

if strcmp(mode,'Batch')
    dataStruct.SpecDif_ind_all = SpecDif_ind_all(t2_idx,:);
    dataStruct.AnalysedDelays  = t2delays(t2_idx);
    return
end

%% Fit and plot the spectral diffusion kinetics
% Define analytical FFCF biexponential decay
    % FFCF:
    %   C(t) = A + B*exp(-t/Tau)
    FFCF_func   = @(P,t) P(1) + P(2)*exp(-t/P(3));

% Define GA fit function
    function fitted_param = runGAfit()
        rng default
%         options = optimoptions('ga','UseParallel', true, 'UseVectorized', false);
        options = optimoptions('ga','FunctionTolerance',1e-15,'PopulationSize',2000);
        fitted_param = ga(@obj_fun,3,[],[],[],[],LB,UB,[]);
        % Nested function to calculate the objectivefunciton
        function ssr = obj_fun(P)
            ssr = sumsqr(FFCF_func(P,t2delays(t2delays>=0))-SpecDif_ind(t2delays>=0));
        end
    end

%% Make the plots and fit the data
% Reload the t2 delays for plotting and fitting
t2delays = dataStruct.t2delays;
% t2delays = t2delays(t2_idx);

switch specdif_options{specdif_typeindx}
    case 'CLS+IvCLS'
        % Fit a bixeponential function to the centerline slope decay, assuming a Kubo lineshape
        LB = [0 0 0];
        UB = [Inf Inf Inf];
        start_param = [1 1 1];
        
        options = optimoptions('lsqcurvefit',...
                    'MaxFunctionEvaluations',1000,...
                    'MaxIterations',1000,...
                    'Algorithm','trust-region-reflective',...  %'levenberg-marquardt' 'trust-region-reflective'
                    'OptimalityTolerance',1e-15,...
                    'FunctionTolerance',1e-15,...
                    'StepTolerance',1e-15,...
                    'Display','off');
        [CLS_param,~,residuals_CLS,~,~,~,jacobian_CLSfit] = lsqcurvefit(FFCF_func,start_param,t2delays(t2_idx),CLS_value(t2_idx),LB,UB,options);
        [IvCLS_param,~,residuals_IvCLS,~,~,~,jacobian_IvCLSfit] = lsqcurvefit(FFCF_func,start_param,t2delays(t2_idx),IvCLS_value(t2_idx),LB,UB,options);
        fittedCLS_CI    = nlparci(CLS_param,residuals_CLS,'jacobian',jacobian_CLSfit);
        fittedCLS_err   = (fittedCLS_CI(:,2) - fittedCLS_CI(:,1))/2;
        fittedIvCLS_CI  = nlparci(IvCLS_param,residuals_IvCLS,'jacobian',jacobian_IvCLSfit);
        fittedIvCLS_err = (fittedIvCLS_CI(:,2) - fittedIvCLS_CI(:,1))/2;
        
        %%%% PLOT CLS
        % Create a new figure with consistent format
        fh = figure();
        fh.Position(3)  = 800;
        fh.Position(4)  = 425;
        fh.Color        = [1 1 1];
        % Define the axes
        ax = axes('Parent',fh);
        ax.FontSize = 14;
        % Plot the spectral diffusion kinetics
        plot(ax,t2delays(t2delays>=0),CLS_value(t2delays>=0),'ob') % 'MarkerFaceColor','b'
        switch fit_method
            case 'None'
                % Do nothing
            otherwise
                % Plot the fit
                FittedSpecDiff   = FFCF_func(CLS_param,linspace(0,max(t2delays)*1.1,500));
                hold on
                plot(ax,linspace(0,max(t2delays),500),FittedSpecDiff,'-b','LineWidth',1)
                hold off
                % Show text with the fit results
                if ~isnan(fittedCLS_err)
                    plot_text = {'Fit of C(t) = A_0 + A_1 exp(-t/\tau) \rm';...
                            ['A_0   = ' num2str(CLS_param(1),'%.3g') ' \pm ' num2str(fittedCLS_err(1),'%.3g')];...
                            ['A_1   = ' num2str(CLS_param(2),'%.3g') ' \pm ' num2str(fittedCLS_err(2),'%.3g')];...
                        ['\bf \tau  = ' num2str(CLS_param(3),'%.3g') ' \pm ' num2str(fittedCLS_err(3),'%.3g') ' ps \rm'];...
                        };
                else
                    plot_text = {'Fit of C(t) = A_0 + A_1 exp(-t/\tau)';...
                        ['A_0 = ' num2str(CLS_param(1),'%.3g')];...
                        ['A_1 = ' num2str(CLS_param(2),'%.3g')];...
                        ['\bf \tau = ' num2str(CLS_param(3),'%.3g') ' ps \rm'];...
                        };
                end
                text(ax,0.65,0.8,plot_text,...
                        'Units','normalized','FontSize',12,'BackgroundColor',textcolor);
        end
        % Make the plot nice
        ylabel(ax,'CLS','FontWeight','bold','FontSize',12);
        xlabel(ax,'Waiting time, t_2 (ps)','FontWeight','bold','FontSize',12);
        ylim(ax,[0,max([max(FittedSpecDiff),max(CLS_value(t2delays>=0))])]);
        % curr_xlim = xlim(ax);
        xlim(ax,[0,max(t2delays)]);
        % Make the axis size consistent
        ax.Units        = 'pixels';
        ax.Position     = [75 70 675 320];
        ax.Units        = 'normalized';
        
        %%% Plot IvCLS
        % Create a new figure with consistent format
        fh = figure();
        fh.Position(3)  = 800;
        fh.Position(4)  = 425;
        fh.Color        = [1 1 1];
        % Define the axes
        ax = axes('Parent',fh);
        ax.FontSize = 14;
        % Plot the spectral diffusion kinetics
        plot(ax,t2delays(t2delays>=0),IvCLS_value(t2delays>=0),'ob') % 'MarkerFaceColor','b'
        switch fit_method
            case 'None'
                % Do nothing
                maxvalue    = max(SpecDif_ind(t2delays>=0));
            otherwise
                % Plot the fit
                FittedSpecDiff   = FFCF_func(IvCLS_param,linspace(0,max(t2delays)*1.1,500));
                hold on
                plot(ax,linspace(0,max(t2delays),500),FittedSpecDiff,'-b','LineWidth',1)
                hold off
                % Show text with the fit results
                if ~isnan(fittedCLS_err)
                    plot_text = {'Fit of C(t) = A_0 + A_1 exp(-t/\tau) \rm';...
                            ['A_0   = ' num2str(IvCLS_param(1),'%.3g') ' \pm ' num2str(fittedIvCLS_err(1),'%.3g')];...
                            ['A_1   = ' num2str(IvCLS_param(2),'%.3g') ' \pm ' num2str(fittedIvCLS_err(2),'%.3g')];...
                        ['\bf \tau  = ' num2str(IvCLS_param(3),'%.3g') ' \pm ' num2str(fittedIvCLS_err(3),'%.3g') ' ps \rm'];...
                        };
                else
                    plot_text = {'Fit of C(t) = A_0 + A_1 exp(-t/\tau)';...
                        ['A_0 = ' num2str(IvCLS_param(1),'%.3g')];...
                        ['A_1 = ' num2str(IvCLS_param(2),'%.3g')];...
                        ['\bf \tau = ' num2str(IvCLS_param(3),'%.3g') ' ps \rm'];...
                        };
                end
                text(ax,0.65,0.8,plot_text,...
                        'Units','normalized','FontSize',12,'BackgroundColor',textcolor);
                maxvalue    = max(FittedSpecDiff);
        end
        % Make the plot nice
        ylabel(ax,'IvCLS','FontWeight','bold','FontSize',12);
        xlabel(ax,'Waiting time, t_2 (ps)','FontWeight','bold','FontSize',12);
        ylim(ax,[0,max([maxvalue,max(IvCLS_value(t2delays>=0))])]);
        % curr_xlim = xlim(ax);
        xlim(ax,[0,max(t2delays)]);
        % Make the axis size consistent
        ax.Units        = 'pixels';
        ax.Position     = [75 70 675 320];
        ax.Units        = 'normalized';
    otherwise
        % Fit a bixeponential function to the centerline slope decay, assuming a Kubo lineshape
        LB = [0 0 0];
        UB = [Inf Inf Inf];
        start_param = [1 1 1];

        switch fit_method
            case 'LSQ'
                options = optimoptions('lsqcurvefit',...
                            'MaxFunctionEvaluations',1000,...
                            'MaxIterations',1000,...
                            'Algorithm','trust-region-reflective',...  %'levenberg-marquardt' 'trust-region-reflective'
                            'OptimalityTolerance',1e-15,...
                            'FunctionTolerance',1e-15,...
                            'StepTolerance',1e-15,...
                            'Display','off');
                [fitted_param,~,residuals,~,~,~,jacobian_fit] = lsqcurvefit(FFCF_func,start_param,t2delays(t2_idx),SpecDif_ind(t2_idx),LB,UB,options);
                fitted_CI   = nlparci(fitted_param,residuals,'jacobian',jacobian_fit);
                fitted_err  = (fitted_CI(:,2) - fitted_CI(:,1))/2;
            case 'GA'
                fitted_param = runGAfit();
                fitted_err  = NaN*ones(1,3);
        end
        % Create a new figure with consistent format
        fh = figure();
        fh.Position(3)  = 800;
        fh.Position(4)  = 425;
        fh.Color        = [1 1 1];

        % Define the axes
        ax = axes('Parent',fh);
        ax.FontSize = 14;

        % Plot the spectral diffusion kinetics
        plot(ax,t2delays(t2delays>=0),SpecDif_ind(t2delays>=0),'ob') % 'MarkerFaceColor','b'
        switch fit_method
            case 'None'
                % Do nothing
                maxvalue    = max(SpecDif_ind(t2delays>=0));
            otherwise
                % Plot the fit
                FittedSpecDiff   = FFCF_func(fitted_param,linspace(0,max(t2delays)*1.1,500));
                hold on
                plot(ax,linspace(0,max(t2delays),500),FittedSpecDiff,'-b','LineWidth',1)
                hold off
                % Show text with the fit results
                if ~isnan(fitted_err)
                    plot_text = {'Fit of C(t) = A_0 + A_1 exp(-t/\tau) \rm';...
                            ['A_0   = ' num2str(fitted_param(1),'%.3g') ' \pm ' num2str(fitted_err(1),'%.3g')];...
                            ['A_1   = ' num2str(fitted_param(2),'%.3g') ' \pm ' num2str(fitted_err(2),'%.3g')];...
                        ['\bf \tau  = ' num2str(fitted_param(3),'%.3g') ' \pm ' num2str(fitted_err(3),'%.3g') ' ps \rm'];...
                        };
                else
                    plot_text = {'Fit of C(t) = A_0 + A_1 exp(-t/\tau)';...
                        ['A_0 = ' num2str(fitted_param(1),'%.3g')];...
                        ['A_1 = ' num2str(fitted_param(2),'%.3g')];...
                        ['\bf \tau = ' num2str(fitted_param(3),'%.3g') ' ps \rm'];...
                        };
                end
                text(ax,0.65,0.8,plot_text,...
                        'Units','normalized','FontSize',12,'BackgroundColor',textcolor);
                maxvalue    = max(FittedSpecDiff);
        end
        
        % Make the plot nice
        ylabel(ax,SpecDif_name,'FontWeight','bold','FontSize',12);
        xlabel(ax,'Waiting time, t_2 (ps)','FontWeight','bold','FontSize',12);
        ylim(ax,[0,max([maxvalue,max(SpecDif_ind(t2delays>=0))])]);
        % curr_xlim = xlim(ax);
        xlim(ax,[0,max(t2delays)]);

        % Make the axis size consistent
        ax.Units        = 'pixels';
        ax.Position     = [75 70 675 320];
        ax.Units        = 'normalized';
end

%% WRITE to dataStruct
% CLS
    dataStruct.CLS_Xdata     = CLS_Xdata;
    dataStruct.CLS_Ydata     = CLS_Ydata;
% IvCLS
    dataStruct.IvCLS_Xdata   = IvCLS_Xdata;
    dataStruct.IvCLS_Ydata   = IvCLS_Ydata;
% NLS
    dataStruct.NLS_Xdata     = NLS_Xdata;
    dataStruct.NLS_Ydata     = NLS_Ydata;

% Show Spectral diffusion control
    app.I2D_ShowSpecDiff.Enable = 'on';
    app.I2D_ShowSpecDiff.Value  = 1;
    dataStruct.SpecDiff         = 1;
    
end
