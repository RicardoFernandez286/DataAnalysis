function SecondaryPlot(handles,where,what)

%% Clear the axes
cla(where,'reset')
hold(where,'off')

%% READ from handles
% Hardcoded values
plotmode            = 'Time';
plot_all            = 0;

% Read the selection on the GUI
m = handles.Population_delay.Value;
k = 1;

% Read common stuff
ProbeAxis           = handles.ProbeAxis;
freq_fit            = handles.freq_fit;
scattering_maxima   = handles.scattering_maxima;
PumpAxis            = handles.PumpAxis{m,k};

% Read phasing data
phased_FFTZPint     = handles.phased_FFTZPint{m,k};
fittedPhase         = handles.fittedPhase{m,k};
ZP_phase            = handles.ZP_phase{m,k};
phasepoints         = handles.phasepoints{m,k};

% Read time-domain data
apodize_function    = handles.apodize_function{m,k};
interferogram       = handles.interferogram{m,k};
signal              = handles.signal{m,k};
% apo_interferogram   = handles.apo_interferogram{m,k};
% apo_signal          = handles.apo_signal{m,k};
t1delays            = handles.t1delays{m,k};
binzero             = handles.binzero{m,k};

%% Plot the selected data type
switch what
    case 'td' % Plot time-domain data for the given pixel
        % Plot the time-domain data
        pixel               = str2double(handles.PixelNumber.String);
        switch plotmode
            case 'Time'
                time = t1delays(:,2);
                xlabel(where,'t_{1} delay (fs)');
            case 'Bins'
                time = t1delays(:,1)+1;
                xlabel(where,'t_{1} delay (bins)');
        end
        if pixel == 0
            y_data = interferogram;
        else
            y_data = signal(:,pixel);
        end
        yyaxis(where,'left')
        plot(where,time,y_data,'-ob','LineWidth',0.5,'MarkerSize',0.5);
        ylabel(where,'Signal (V / mOD)')
        axis(where,'tight')
        hold(where,'on');
        % Plot the apodization function and bin zero
        yyaxis(where,'right')
        plot(where,time,apodize_function,'-r','LineWidth',2);
        plot(where,time(binzero),apodize_function(binzero),'or','LineWidth',2);
        ylabel(where,'Apo. Func.');
        ylim(where,[-1,1]);
        hold off
        xlim(where,[min(time) max(time)]);
        zl = refline(where,0,0);
        zl.Color = [0.5 0.5 0.5];
    case 'ph' % Plot phasing data
        %%% Plot the pump spectrum (real part of FFT[pyro])
        yyaxis(where,'left')
        p1 = plot(where,PumpAxis,abs(phased_FFTZPint),'-b','LineWidth',2,'DisplayName','|\mathcal{F}(Int)|');
        hold(where,'on')
        if plot_all == 1
            % Plot also the real and imaginary parts of the pump spectrum
            p2 = plot(where,PumpAxis,real(phased_FFTZPint),'-c','DisplayName','Re[\mathcal{F}(Int)]');
            p3 = plot(where,PumpAxis,imag(phased_FFTZPint),'-','Color',[0.9290, 0.6940, 0.1250],'DisplayName','Im[\mathcal{F}(Int)]');
        end
        ylabel(where,'Intensity (a.u.)');
        hline = refline(where,0,0);
        hline.Color = [0.5 0.5 0.5];
        hline.LineWidth = 0.1;
        %%% Plot the phase and the fit
        yyaxis(where,'right')
        hold(where,'on')
        % Fitted points
        p4 = plot(where,PumpAxis(phasepoints),fittedPhase(phasepoints),'dg','LineWidth',0.1,'MarkerSize',6,'DisplayName','Fitted phase');
        p5 = plot(where,PumpAxis,fittedPhase,'g','LineWidth',3,'DisplayName','Fitted phase');
        % Original phase
        p6 = plot(where,PumpAxis,ZP_phase,'-r','LineWidth',1,'DisplayName','Calculated phase');
        ylabel(where,'Phase (rad)')
        % Adjust the limits (Y axis - phase)
        ylim(where,[-pi,pi])
        yticks(where,[-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi ])
        yticklabels(where,{'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi'})
        % Adjust the limits (X axis)
        xlim(where,[min(ProbeAxis)*0.95 max(ProbeAxis)*1.05]);
        xlabel(where,'Pump wavenumber (cm^{-1})');
        if plot_all == 1
            % Show legend
            legend(where,[p1 p2 p3],{"$\mid\hat{\mathcal{F}}$(Int)$\mid$" "Re[$\hat{\mathcal{F}}$(Int)]" "Im[$\hat{\mathcal{F}}$(Int)]"},'Interpreter','latex');
            legend(where,'Boxoff')
            legend(where,'Location','northwest')
        end
    case 'pc' % Plot probe calibration (data+fit)
        plot(where,1:length(ProbeAxis),scattering_maxima,'-xr','LineWidth',2);
        xlabel(where,'Pixel number');
        ylabel(where,'Fitted freq. (cm^{-1})');
        hold(where,'on')
        plot(where,1:length(ProbeAxis),freq_fit,'-k','LineWidth',1);
        hold(where,'off')
        xlim(where,[1,length(ProbeAxis)]);
        ylim(where,[min(freq_fit),max(freq_fit)]);
        % Show legend
        legend(where,{'Scattering maxima' 'Fit'})
        legend(where,'Boxoff')
        legend(where,'Location','northwest')
end