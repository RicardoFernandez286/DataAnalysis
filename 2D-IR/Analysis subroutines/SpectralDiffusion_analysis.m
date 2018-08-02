function  handles = SpectralDiffusion_analysis(handles)

% Description: This function handles the spectral diffusion analysis of a
% series of 2D-IR data as a function of the waiting time
% 
% Usage: handles = SpectralDiffusion_analysis(handles)
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
% Ricardo Fernández-Terán / 17.07.2018 / v1.0a

%% READ from handles
% Get information from the GUI
plot_pumpdirection  = char(handles.plot_pumpdirection.String{handles.plot_pumpdirection.Value});
ProbeAxis           = handles.ProbeAxis;
PumpAxis            = handles.PumpAxis;
Ndelays             = handles.Ndelays;
t2delays            = handles.t2delays;
FFT_ZPsig           = handles.FFT_ZPsig;
PROC_2D_DATA        = handles.PROC_2D_DATA;
Npixels             = length(ProbeAxis);
Nbins               = length(PumpAxis{1,1});
Axes                = handles.MainAxes;
popdelay            = handles.Population_delay.Value;

% Hardcoded settings
Interp_method       = 'InterpFT'; % 2D FFT, InterpFT, Mesh 
Interp_order        = 4; % Multiplies Npixels by a factor to interpolate in the probe dimension
% Fit_type            = 'Quadratic'; % 'Quadratic' 'None' - Sets the type of fit to calculate the minima along the slices
N_pointsParabola    = 10; % No. of X points to fit a parabola on each side of the peak: [-x (peak) +x]
intensity_threshold = 60/100; % Intensity threshold for the linear fit to get the CLS or IvCLS
textcolor           = 'none'; % 'none' or RGB color
fit_method          = 'LSQ'; % 'GA' or 'LSQ' for Genetic algorithm or Least-squares

% Suppress warnings
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');

% Initialize variables (ensures only the current analysis is saved)
CLS_Xdata     = NaN;
CLS_Ydata     = NaN;
IvCLS_Xdata   = NaN;
IvCLS_Ydata   = NaN;
NLS_Xdata     = NaN;
NLS_Ydata     = NaN;

%% Ask the user which method to use
specdif_options     = {'CLS','IvCLS','NLS','CLS+IvCLS','Gaussian Fit'};
[specdif_typeindx,doAnalysis] = listdlg('ListString',specdif_options,'OKstring','Continue','SelectionMode','single','ListSize',[150,100],'PromptString','Select analysis type:');

if doAnalysis == 0
    return
end

%% Interpolate the data along the probe axis (interpolation along w1 can be performed by apodisation)
switch Interp_method
    case '2D FFT'
        for m=1:Ndelays
            Interp_Data{m,1}    = real(fft(ifft(PROC_2D_DATA{m,1},2),Npixels.*Interp_order,2));
        end
    case 'InterpFT'
        for m=1:Ndelays
            Interp_Data{m,1}    = interpft(PROC_2D_DATA{m,1},Npixels.*Interp_order,2);
        end
    case 'Mesh'
        % Do something
end

%% Get the spectral diffusion dynamics
% Interpolate the probe axis
pixels              = 1:1:Npixels;
Newpixels           = linspace(1,32,Npixels.*Interp_order);
Interp_ProbeAxis    = interp1(pixels,ProbeAxis,Newpixels);

switch handles.InteractiveModeTick.Value
    case 1
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
    case 0
        RegionPositions_cell    = inputdlg({'Enter the pump regions (Format: min1 max1; min2 max2; ...): ','Enter the probe regions (Format: min1 max1; min2 max2; ...):'},'Define regions to integrate', [1 80]);
        pump_ranges             = str2num(RegionPositions_cell{1});
        probe_ranges            = str2num(RegionPositions_cell{2});
        L = size(pump_ranges,1);
        M = size(probe_ranges,1);
        if L ~= M
            warndlg('Number of pump and probe ranges is different!')
            return
        end
end

% Parse the pump ranges into indexes
pump_idxrange  = findClosestId2Val(PumpAxis{1,1},pump_ranges);
probe_idxrange = findClosestId2Val(Interp_ProbeAxis,probe_ranges);

% Get the data according to the selected method
switch specdif_options{specdif_typeindx}
    case 'CLS'
        % Minima in w3 for each value of w1
        probe_indexes       = probe_idxrange(1):1:probe_idxrange(end);
        pump_indexes        = pump_idxrange(1):1:pump_idxrange(end);
        Npeaks              = length(pump_indexes);
        min_values          = zeros(Npeaks,Ndelays);
        SpecDif_ind         = zeros(Ndelays,1);
        for m=1:Ndelays
            % Get the data around the minimum
            [peak_val,peak_pos] = min(Interp_Data{m,1}(pump_indexes,probe_indexes),[],2);
            peak_pos            = peak_pos + probe_idxrange(1);
            % Fit a quadratic polynomial within 2 points around the minimum 
            %   and find the analytical vertex of the parabola
            for i=1:Npeaks
                probe_segment   = Interp_ProbeAxis(peak_pos(i)-N_pointsParabola:peak_pos(i)+N_pointsParabola);
                data_segment    = Interp_Data{m,1}(pump_indexes(i),peak_pos(i)-N_pointsParabola:peak_pos(i)+N_pointsParabola);
                coeff           = polyfit(probe_segment,data_segment,2);
                min_values(i,m) = -coeff(2)/(2*coeff(1));
            end
            % Fit a line around the minimum according to the intensity cutoff
            peak_val            = peak_val/max(abs(peak_val));
            above_threshold     = abs(peak_val) >= intensity_threshold;
            pump_cut            = pump_indexes(above_threshold);
            mdl                 = fit(PumpAxis{1,1}(pump_cut),min_values(above_threshold,m),'poly1','Robust','Bisquare');
            CLS_coeff           = coeffvalues(mdl);
            SpecDif_ind(m)      = CLS_coeff(1);
        end
        SpecDif_name = 'CLS';
%         hold on
%         plot(PumpAxis{1,1}(pump_cut),min_values(above_threshold,popdelay),'Color','w','LineWidth',2)
%         hold off
        % Save the spectral diffusion data
        CLS_Xdata     = PumpAxis{1,1}(pump_cut);
        CLS_Ydata     = min_values(above_threshold,:);
    case 'IvCLS'
        % Minima in w1 for each value of w3
        probe_indexes       = probe_idxrange(1):1:probe_idxrange(end);
        pump_indexes        = pump_idxrange(1):1:pump_idxrange(end);
        Npeaks              = length(probe_indexes);
        min_values          = zeros(Npeaks,Ndelays);
        SpecDif_ind         = zeros(Ndelays,1);
        for m=1:Ndelays
            % Get the data around the minimum
            [peak_val,peak_pos] = min(Interp_Data{m,1}(pump_indexes,probe_indexes),[],1);
            peak_pos            = peak_pos + pump_idxrange(1);
            % Fit a quadratic polynomial within 2 points around the minimum 
            %   and find the analytical vertex of the parabola
            for i=1:Npeaks
                pump_segment    = PumpAxis{1,1}(peak_pos(i)-N_pointsParabola:peak_pos(i)+N_pointsParabola);
                data_segment    = Interp_Data{m,1}(peak_pos(i)-N_pointsParabola:peak_pos(i)+N_pointsParabola,probe_indexes(i));
                coeff           = polyfit(pump_segment,data_segment,2);
                min_values(i,m) = -coeff(2)/(2*coeff(1));
            end
            % Fit a line around the minimum according to the intensity cutoff
            peak_val            = peak_val/max(abs(peak_val));
            above_threshold     = abs(peak_val) >= intensity_threshold;
            probe_cut           = probe_indexes(above_threshold);
            mdl                 = fit(min_values(above_threshold,m),Interp_ProbeAxis(probe_cut)','poly1','Robust','Bisquare');
            IvCLS_coeff         = coeffvalues(mdl);
            SpecDif_ind(m)      = 1./IvCLS_coeff(1);
        end
        SpecDif_name    = 'IvCLS';
        % Save the spectral diffusion data
        IvCLS_Xdata     = min_values(above_threshold,:);
        IvCLS_Ydata     = Interp_ProbeAxis(probe_cut);
    case 'NLS'
        
    case 'CLS+IvCLS'
        % IvCLS
        % Minima in w1 for each value of w3
        probe_indexes       = probe_idxrange(1):1:probe_idxrange(end);
        pump_indexes        = pump_idxrange(1):1:pump_idxrange(end);
        Npeaks              = length(probe_indexes);
        min_values          = zeros(Npeaks,Ndelays);
        IvCLS_value         = zeros(Ndelays,1);
        for m=1:Ndelays
            % Get the data around the minimum
            [peak_val,peak_pos] = min(Interp_Data{m,1}(pump_indexes,probe_indexes),[],1);
            peak_pos            = peak_pos + pump_idxrange(1);
            % Fit a quadratic polynomial within 2 points around the minimum 
            %   and find the analytical vertex of the parabola
            for i=1:Npeaks
                pump_segment    = PumpAxis{1,1}(peak_pos(i)-N_pointsParabola:peak_pos(i)+N_pointsParabola);
                data_segment    = Interp_Data{m,1}(peak_pos(i)-N_pointsParabola:peak_pos(i)+N_pointsParabola,probe_indexes(i));
                coeff           = polyfit(pump_segment,data_segment,2);
                min_values(i,m) = -coeff(2)/(2*coeff(1));
            end
            % Fit a line around the minimum according to the intensity cutoff
            peak_val            = peak_val/max(abs(peak_val));
            above_threshold     = abs(peak_val) >= intensity_threshold;
            probe_cut           = probe_indexes(above_threshold);
            mdl                 = fit(min_values(above_threshold,m),Interp_ProbeAxis(probe_cut)','poly1','Robust','Bisquare');
            IvCLS_coeff         = coeffvalues(mdl);
            IvCLS_value(m)      = 1./IvCLS_coeff(1);
        end
        SpecDif_name    = 'IvCLS';
        % Save the spectral diffusion data
        IvCLS_Xdata     = min_values(above_threshold,:);
        IvCLS_Ydata     = Interp_ProbeAxis(probe_cut);
        
        % CLS
        % Minima in w3 for each value of w1
        probe_indexes       = probe_idxrange(1):1:probe_idxrange(end);
        pump_indexes        = pump_idxrange(1):1:pump_idxrange(end);
        Npeaks              = length(pump_indexes);
        min_values          = zeros(Npeaks,Ndelays);
        CLS_value           = zeros(Ndelays,1);
        for m=1:Ndelays
            % Get the data around the minimum
            [peak_val,peak_pos] = min(Interp_Data{m,1}(pump_indexes,probe_indexes),[],2);
            peak_pos            = peak_pos + probe_idxrange(1);
            % Fit a quadratic polynomial within 2 points around the minimum 
            %   and find the analytical vertex of the parabola
            for i=1:Npeaks
                probe_segment   = Interp_ProbeAxis(peak_pos(i)-N_pointsParabola:peak_pos(i)+N_pointsParabola);
                data_segment    = Interp_Data{m,1}(pump_indexes(i),peak_pos(i)-N_pointsParabola:peak_pos(i)+N_pointsParabola);
                coeff           = polyfit(probe_segment,data_segment,2);
                min_values(i,m) = -coeff(2)/(2*coeff(1));
            end
            % Fit a line around the minimum according to the intensity cutoff
            peak_val            = peak_val/max(abs(peak_val));
            above_threshold     = abs(peak_val) >= intensity_threshold;
            pump_cut            = pump_indexes(above_threshold);
            mdl                 = fit(PumpAxis{1,1}(pump_cut),min_values(above_threshold,m),'poly1','Robust','Bisquare');
            CLS_coeff           = coeffvalues(mdl);
            CLS_value(m)        = CLS_coeff(1);
        end
        SpecDif_name = 'CLS';
        % Save the spectral diffusion data
        CLS_Xdata     = PumpAxis{1,1}(pump_cut);
        CLS_Ydata     = min_values(above_threshold,:);
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
                    'StepTolerance',1e-15);
        [CLS_param,~,residuals_CLS,~,~,~,jacobian_CLSfit] = lsqcurvefit(FFCF_func,start_param,t2delays(t2delays>=0),CLS_value(t2delays>=0),LB,UB,options);
        [IvCLS_param,~,residuals_IvCLS,~,~,~,jacobian_IvCLSfit] = lsqcurvefit(FFCF_func,start_param,t2delays(t2delays>=0),IvCLS_value(t2delays>=0),LB,UB,options);
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
        ylim(ax,[0,max(FittedSpecDiff)]);
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
        ylim(ax,[0,maxvalue]);
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
                            'StepTolerance',1e-15);
                [fitted_param,~,residuals,~,~,~,jacobian_fit] = lsqcurvefit(FFCF_func,start_param,t2delays(t2delays>=0),SpecDif_ind(t2delays>=0),LB,UB,options);
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
        ylim(ax,[0,maxvalue]);
        % curr_xlim = xlim(ax);
        xlim(ax,[0,max(t2delays)]);

        % Make the axis size consistent
        ax.Units        = 'pixels';
        ax.Position     = [75 70 675 320];
        ax.Units        = 'normalized';
end

%% WRITE to handles
% CLS
    handles.CLS_Xdata     = CLS_Xdata;
    handles.CLS_Ydata     = CLS_Ydata;
% IvCLS
    handles.IvCLS_Xdata   = IvCLS_Xdata;
    handles.IvCLS_Ydata   = IvCLS_Ydata;
% NLS
    handles.NLS_Xdata     = NLS_Xdata;
    handles.NLS_Ydata     = NLS_Ydata;

% Show Spectral diffusion control
    handles.ShowSpecDiff.Visible = 'on';
    handles.SpecDiff      = 1;
    
end
