function [cm,lam] = CalibrateWLAxis(CAL_data)
% Subroutine to calibrate the wavelength axis of a spectrometer, using a centering and scaling approach.
% Will output the fitted cm-1 and nm axes
% The number of columns of cm and lam = Ndet (given by the size of the AbsSolvent cell in the CAL_data structure).
% Ricardo Fernandez-Teran / 24.06.2022 / v3.0c

doFit = 1;
%% Read from input structure
CalType     = CAL_data.CalType;

if CalType == 1
    Ndet = 2;
else
    Ndet = 1;
end

MeasY_all   = CAL_data.AbsSolvent;

for i=1:Ndet
    MeasY   = MeasY_all{i};
    MeasX   = 1:length(MeasY);

    RefX    = CAL_data.CalX;
    RefY    = CAL_data.CalY;

    % Ensure everything is stored as column vectors
    MeasX   = reshape(MeasX,[],1);
    MeasY   = reshape(MeasY,[],1);
    RefX    = reshape(RefX,[],1);
    RefY    = reshape(RefY,[],1);

    % If we're doing IR calibration, convert to nm since the grating is linear in nm
    % The X axis of the reference spectrum is converted to nm, and the Y experimental data is flipped around
    switch CalType
        case 1 % UoS 2DIR/TRIR
            Gratings    = CAL_data.Gratings;
            CWL         = CAL_data.CWL;
            RefX        = 1e7./RefX;
            MeasY       = flipud(MeasY);
        case 4 % UZH Lab 2
            Gratings    = CAL_data.Gratings;
            CWL         = CAL_data.CWL;
            RefX        = 1e7./RefX;
            MeasY       = flipud(MeasY);
    end

    %% Cut off relevant part of reference spectrum -- set initial guesses for fit parameters
    switch CalType
        case 1 % UoS 2DIR/TRIR
            switch Gratings(i)
                case 0 % 120 l/mm
                    minRefX = CWL(i) - 200;
                    maxRefX = CWL(i) + 200;
                    ppnm    = 0.153;             % pixels per nm = 1/resolution;
                case 1 % 100 l/mm
                    minRefX = CWL(i) - 350;
                    maxRefX = CWL(i) + 350;
                    ppnm    = 0.137;             % pixels per nm = 1/resolution;
                case 2 % 50 l/mm
                    minRefX = CWL(i) - 350;
                    maxRefX = CWL(i) + 550;
                    ppnm    = 0.068;             % pixels per nm = 1/resolution;
            end
            wp1         = CWL(i) - 48/ppnm; % lowest wavelength (nm) = CWL - Npix/2*resolution
            cutID       = sort(findClosestId2Val(RefX,[minRefX maxRefX]));
            RefX_cut    = RefX(cutID(1):cutID(2));
            RefY_cut    = RefY(cutID(1):cutID(2));
        
        case 2    % UniGE TA 
            minRefX  = 345;
            maxRefX  = 690;

            cutID    = sort(findClosestId2Val(RefX,[minRefX maxRefX]));
            [RefX_cut,srt_id] = sort(RefX(cutID(1):cutID(2)));
            RefY_cut = RefY(cutID(1):cutID(2));
            RefY_cut = RefY_cut(srt_id);

            RefY     = RefY./max(RefY_cut);
            RefY_cut = RefY_cut./max(RefY_cut);
            
        case 3     % UniGE nsTA 
            minRefX  = 370;
            maxRefX  = 690;

            cutID    = sort(findClosestId2Val(RefX,[minRefX maxRefX]));
            [RefX_cut,srt_id] = sort(RefX(cutID(1):cutID(2)));
            RefY_cut = RefY(cutID(1):cutID(2));
            RefY_cut = RefY_cut(srt_id);

            RefY     = RefY./max(RefY_cut);
            RefY_cut = RefY_cut./max(RefY_cut);
        case 4     % UZH Lab 2
            switch Gratings(i) 
                case 100 % 100 l/mm
                    minRefX = CWL(i) - 350;
                    maxRefX = CWL(i) + 350;
                    ppnm    = 0.137;             % pixels per nm = 1/resolution;
                case 150 % 150 l/mm
                    minRefX = CWL(i) - 150;
                    maxRefX = CWL(i) + 150;
                    ppnm    = 0.068;             % pixels per nm = 1/resolution;
                case 300
                    minRefX = CWL(i) - 100;
                    maxRefX = CWL(i) + 100;
                    ppnm    = 0.034;             % pixels per nm = 1/resolution;
            end
            wp1         = CWL(i) - 32/ppnm;      % lowest wavelength (nm) = CWL - Npix/2*resolution
            cutID       = sort(findClosestId2Val(RefX,[minRefX maxRefX]));
            RefX_cut    = RefX(cutID(1):cutID(2));
            RefY_cut    = RefY(cutID(1):cutID(2));
        case 5  % UniGE NIR-TA
            minRefX  = 800;
            maxRefX  = 1450;

            cutID    = sort(findClosestId2Val(RefX,[minRefX maxRefX]));
            [RefX_cut,srt_id] = sort(RefX(cutID(1):cutID(2)));
            RefY_cut = RefY(cutID(1):cutID(2));
            RefY_cut = RefY_cut(srt_id);

            RefY     = RefY./max(RefY_cut);
            RefY_cut = RefY_cut./max(RefY_cut);
            
            MeasX    = (1:512)';
            MeasY    = MeasY(MeasX);
    end

    %% Define Fit Functions
    ref_fun     = griddedInterpolant(RefX_cut,RefY_cut);
    meas_fun    = griddedInterpolant(MeasX,MeasY);
    switch CalType
        case 5 % UniGE NIR-TA
            % Refractive index -- Sellmeier equation (lambda in µm)
            n_l     = @(L,b,c) sqrt(1+sum(b.*(L./1000).^2./((L./1000).^2-c),2));
            
            % Prism properties
            A       = deg2rad(60.84);   % apex angle of the prism
            WL_r    = 1100;             % Reference wavelength (nm)
            
            % Sellmeier coefficients for SF10 (lambda in µm)
            SF10.B  = [1.62153902    0.256287842     1.64447552];
            SF10.C  = [0.0122241457  0.0595736775    147.468793];

            % Refractive index of SF10
            nSF10   = @(l) n_l(l,SF10.B,SF10.C);
            
            % Output angle from the prism
            th_out  = @(l,th_in) asin(nSF10(l).*sin(A - asin(deg2rad(th_in)./nSF10(l))));
            
            % Parameters: 1=focal distance(mm) 2=Xshift 3=input angle  4=Yscale 5=Yshift            
            pix_fun     = @(p) p(1).*1000./25.*sin(th_out(RefX_cut,p(3)) - th_out(WL_r,p(3))) + p(2);
            plot_fun    = @(p) (meas_fun(pix_fun(p))+p(5)).*p(4);
            fit_fun     = @(p) (plot_fun(p) - ref_fun(RefX_cut)).*1e4;
        otherwise
            % Parameters: 1=centralWL 2=nm/pix 3=Yscale 4=Yshift            
            plot_fun    = @(p) meas_fun((RefX_cut-p(1)).*p(2)).*p(3)+p(4)+1e-6.*p(5).*(RefX_cut-p(1))+1e-6.*p(6).*(RefX_cut-p(1)).^2;
            fit_fun     = @(p) (plot_fun(p) - ref_fun(RefX_cut)).*1e4;
    end
    %% Calibration Fit

    switch CalType
        case 1 % UoS TRIR/2DIR
            p0 = [wp1  ppnm  1.5   -0.1  eps  eps  ];
            LB = [1e3  1e-4  0.1  -1    -Inf -Inf ];
            UB = [1e4  0.25  10   1     Inf  Inf  ];
            ftol= 5e-7;
            stol= 5e-7;
        case 2 % UniGE TA
            p0 = [310 1.05  0.6 -0.1 1 1];
            LB = [150 0.5   0.1 -Inf -Inf -Inf];
            UB = [400 2     10   Inf  Inf Inf ];
            ftol= 5e-7;
            stol= 5e-7;
        case 3 % UniGE nsTA
            p0 = [330 1.33  0.58 -0.1 1     1];
            LB = [100 0.1   0.1  -1 -Inf -Inf];
            UB = [400 2     10   1  Inf  Inf ];
            ftol= 1e-6;
            stol= 1e-6;
        case 4 % UZH Lab 2
            p0 = [wp1  ppnm  1    -0.1  eps eps];
            LB = [1e3  1e-4  0.1  -1    -1  -1 ];
            UB = [1e4  0.25  10   1     1   1  ];
            ftol= 5e-7;
            stol= 5e-7;
        case 5 % UniGE NIR-TA
            p0 = [248 175 60  2   -0.5  ];
            LB = [100 100 30  0.1 -1    ];
            UB = [370 240 90  5   +1    ];
            ftol= 1e-12;
            stol= 1e-12;
    end
    
    if doFit == 1
        mb = msgbox('Fit in progress!','Fitting...','help');
        
        options     = optimoptions(@lsqnonlin,...
                    'FunctionTolerance',ftol,...
                    'MaxIterations',5e3,...
                    'MaxFunctionEvaluations',5e4,...
                    'steptolerance',stol,...
                    'display','off',...
                    'typicalX',p0,...
                    'plotfcn','optimplotresnorm_log');
    
        Pfit(i,:) = lsqnonlin(fit_fun,p0,LB,UB,options);
        delete(mb);
    else
        Pfit(i,:)=p0;
    end

    switch CalType
        case {1,4} % IR Setups UoS and UZH Lab 2
            lam(:,i)    = flipud(MeasX/Pfit(i,2) + Pfit(i,1));
        case 5
            % Need to numerically invert the pixel(lambda) function to get lambda(pixel)
            wl_map      = 1e7./(5000:0.1:20000)';
            pix_fun_map = @(p) p(1).*1000./25.*sin(th_out(wl_map,p(3)) - th_out(WL_r,p(3))) + p(2);
            idxAll      = (1:length(MeasY_all{1}))';
            lam_map     = pix_fun_map(Pfit(i,:));
            lam(:,i)    = wl_map(knnsearch(lam_map,idxAll));
        otherwise 
            lam(:,i)    = MeasX/Pfit(i,2) + Pfit(i,1); 
    end
    cm(:,i)         = 1e7./lam(:,i);

    %% Plot Fit Results
    fh(i) = figure(i); %#ok<*AGROW> 
    fh(i).Color = 'w';
    fh(i).Position(3:4) = [600 350];
    fh(i).Name = ['Detector No.' num2str(i) ' Calibration / Exp(Red) vs Ref(Black)'];
    clf(fh(i));
    
    ax = axes('parent',fh(i));

    if i>1
        fh(i).Position(2) = fh(i-1).Position(2) - fh(i-1).Position(4) - 100;
    end

    switch CalType
        case {1,4} % IR setups UoS and UZH Lab 2
            plot(ax,1e7./RefX,RefY,'k')
            hold(ax,'on');
            plot(ax,1e7./RefX_cut,plot_fun(Pfit(i,:)),'r','LineWidth',2)
            hold(ax,'off');

            axis(ax,'tight');
            xlim(ax,[min(1e7./RefX_cut)-20 max(1e7./RefX_cut)+20]);
            xlabel(ax,'Real Wavenumbers (cm^{-1})','FontWeight','bold');
            ylabel(ax,'Norm. Absorbance','FontWeight','bold');
            title(ax,['Detector ' num2str(i)],'FontWeight','bold');
            ax.FontSize = 16;
        
        case {2,3,5}
            plot(ax,RefX,RefY,'k')
            hold(ax,'on');
            plot(ax,RefX_cut,plot_fun(Pfit(i,:)),'r','LineWidth',2)
            hold(ax,'off');

            axis(ax,'tight');
            xlim(ax,[min(RefX_cut)-10 max(RefX_cut)+10]);
            xlabel(ax,'Real Wavelength (nm)','FontWeight','bold');
            ylabel(ax,'Norm. Absorbance','FontWeight','bold');
            title(ax,['Detector ' num2str(i)],'FontWeight','bold');
            ax.FontSize = 16;
    end
end