function SecondaryPlot(dataStruct,where,what,varargin)

if dataStruct.isSimulation == 1
    return
end

%% Clear the axes
cla(where,'reset')
hold(where,'off')

%% READ from handles
% Hardcoded values
plotmode            = 'Time';
plot_all            = 0;

% Read the selection on the GUI
m = dataStruct.plotOptions.popdelay;
k = 1;

% Read common stuff
ProbeAxis           = dataStruct.ProbeAxis;
freq_fit            = dataStruct.freq_fit;
scattering_maxima   = dataStruct.scattering_maxima;
PumpAxis            = dataStruct.PumpAxis{m,k};

% Read phasing data
phased_FFTZPint     = dataStruct.phased_FFTZPint{m,k};
fittedPhase         = dataStruct.fittedPhase{m,k};
ZP_phase            = dataStruct.ZP_phase{m,k};
phasepoints         = dataStruct.phasepoints{m,k};

% Read time-domain data
apodize_function    = dataStruct.apodize_function{m,k};
interferogram       = dataStruct.interferogram{m,k};
signal              = dataStruct.signal{m,k};
if ~isempty(varargin)
    pixel           = varargin{1};
    pixline         = varargin{2};
    pixline_axis    = varargin{3};
end

% apo_interferogram   = dataStruct.apo_interferogram{m,k};
% apo_signal          = dataStruct.apo_signal{m,k};
t1delays            = dataStruct.t1delays{m,k};
binzero             = dataStruct.binzero{m,k};

%% Plot the selected data type
switch what
    case 'td' % Plot time-domain data for the given pixel
        % Plot the time-domain data
        switch plotmode
            case 'Time'
                time = t1delays(:,2);
                xlabel(where,'t_{1} delay (fs)','FontWeight','bold');
            case 'Bins'
                time = t1delays(:,1)+1;
                xlabel(where,'t_{1} delay (bins)','FontWeight','bold');
        end
        if pixel == 0
            y_data = interferogram;
            label  = ' V';
        else
            y_data = signal(:,pixel);
            label  = ' mOD';
        end
        yyaxis(where,'left')
        plot(where,time,y_data,'-ob','LineWidth',0.5,'MarkerSize',0.5);
        ylabel(where,'Signal (V / mOD)','FontWeight','bold')
        axis(where,'tight')
        hold(where,'on');
        % Plot the apodization function and bin zero
        yyaxis(where,'right')
        plot(where,time,apodize_function,'-r','LineWidth',2);
        plot(where,time(binzero),apodize_function(binzero),'or','LineWidth',2);
        ylabel(where,'Apo. Func.','FontWeight','bold');
        ylim(where,[-1,1]);
        hold(where,'off')
        xlim(where,[min(time) max(time)]);
        zl = yline(where,0);
        zl.Color = [0.5 0.5 0.5];
        % Show info about maximum/minimum values
        maxvalue    = num2str(max(y_data),'%.3g');
        minvalue    = num2str(min(y_data),'%.3g');
        text(where,0.5,0.15,{['Max = ' maxvalue label];['Min = ' minvalue label]},'Units','normalized','FontSize',14);
%         % Show interferometer contrast (wrong formula! needs to be checked!)
%         if pixel == 0
%             contrast    = num2str(100.*(max(y_data)-min(y_data)-mean(y_data)-mean(abs([max(y_data),min(y_data)])))./(max(y_data)-min(y_data)+mean(y_data)+mean(abs([max(y_data),min(y_data)]))),'%.3g');
%             text(where,0.635,0.9,['Contrast: ' contrast '%'],'Units','normalized','FontSize',10);
%         end
        yyaxis(where,'left')
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
        ylabel(where,'Intensity (a.u.)','FontWeight','bold');
        hline = yline(where,0);
        hline.Color = [0.5 0.5 0.5];
        hline.LineWidth = 0.1;
        %%% Plot the phase and the fit
        yyaxis(where,'right')
        hold(where,'on')
        % Fitted points
        p4 = plot(where,PumpAxis(phasepoints),wrapToPi(fittedPhase(phasepoints)),'dg','LineWidth',0.1,'MarkerSize',6,'DisplayName','Fitted phase');
        p5 = plot(where,PumpAxis,wrapToPi(fittedPhase),'g','LineWidth',3,'DisplayName','Fitted phase');
        % Original phase
        p6 = plot(where,PumpAxis,wrapToPi(ZP_phase),'-r','LineWidth',1,'DisplayName','Calculated phase');
        ylabel(where,'Phase (rad)','FontWeight','bold')
        % Adjust the limits (Y axis - phase)
        ylim(where,[-pi,pi])
        yticks(where,[-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi ])
        yticklabels(where,{'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi'})
        % Adjust the limits (X axis)
        xlim(where,[min(ProbeAxis)*0.95 max(ProbeAxis)*1.05]);
        xlabel(where,'Pump wavenumber (cm^{-1})','FontWeight','bold');
        if plot_all == 1
            % Show legend
            legend(where,[p1 p2 p3],["$\mid\hat{\mathcal{F}}$(Int)$\mid$" "Re[$\hat{\mathcal{F}}$(Int)]" "Im[$\hat{\mathcal{F}}$(Int)]"],'Interpreter','latex','FontSize',14);
            legend(where,'boxoff')
            legend(where,'Location','best')
        end
    case 'pc' % Plot probe calibration (data+fit)
        plot(where,1:length(ProbeAxis),scattering_maxima,'-xr','LineWidth',2);
        xlabel(where,'Pixel number','FontWeight','bold');
        ylabel(where,'Fitted freq. (cm^{-1})','FontWeight','bold');
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


if pixline == 1 && pixel > 0
    if strcmp(class(pixline_axis.Children(1)),'matlab.graphics.chart.decoration.ConstantLine')
        pixline_axis.Children(1).Value = ProbeAxis(pixel);
    else
        yl          = yline(pixline_axis,ProbeAxis(pixel));
        yl.Color    = 0.5*[1 1 1];
        yl.LineWidth= 3;
    end
elseif pixline == 0
    if strcmp(class(pixline_axis.Children(1)),'matlab.graphics.chart.decoration.ConstantLine')
        delete(pixline_axis.Children(1));
    end
end
        
where.Box = 'on';