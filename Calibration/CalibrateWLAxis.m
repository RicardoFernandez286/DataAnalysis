function [cm,lam] = CalibrateWLAxis(CAL_data,FitLimits,doFit)
% Subroutine to calibrate the wavelength axis of a spectrometer, using a centering and scaling approach.
% Will output the fitted cm-1 and nm axes
% The number of columns of cm and lam = Ndet (given by the size of the AbsSolvent cell in the CAL_data structure).
% Ricardo Fernandez-Teran / 08.02.2023 / v3.2d

%% Read from input structure
CalType     = CAL_data.CalType;

switch CalType
    case {1,6,7}; Ndet = 2;
    otherwise;    Ndet = 1;
end

MeasY_all   = CAL_data.AbsSolvent;

for i=1:Ndet
    MeasY   = MeasY_all{i};
    Npix(i) = length(MeasY);
    MeasX   = 1:Npix(i);
    
    % Reference (standard) spectrum
    RefX    = CAL_data.CalX;
    RefY    = CAL_data.CalY;

    % Ensure everything is stored as column vectors
    MeasX   = reshape(MeasX,[],1);
    MeasY   = reshape(MeasY,[],1);
    RefX    = reshape(RefX,[],1);
    RefY    = reshape(RefY,[],1);
    
    switch CalType
        % If we're doing IR calibration, start from high to low wavenumbers
        % (i.e. low to high wavelength!!)
        case {1,4,6,7}
            [RefX,idRefX]    = sort(RefX,'descend');
            RefY             = RefY(idRefX);
    end

    % If we're doing IR calibration, convert standard spectral axis to nm since the grating is linear in nm
    % The X axis of the reference [standard] spectrum is converted to nm, and the Y experimental data is flipped around
    switch CalType
        case {1,4} % {UoS 2DIR/TRIR, UZH Lab 1,2,4, UniGE TRIR new}
            Gratings    = CAL_data.Gratings;
            CWL         = CAL_data.CWL;
            RefX        = 1e7./RefX;
            MeasY       = flipud(MeasY);
        case {2,3} % UniGE 
            CWL         = FitLimits{3,1};
            maxWL       = FitLimits{1,1};
            minWL       = FitLimits{2,1};
            ppnmT       = FitLimits{4,1};
        case {6,7} % {RAL Abs, Ral Int}
            CWL         = 1e7./[FitLimits{3,1:2}];
            maxWL       = 1e7./[FitLimits{1,1:2}]; % max WL = min cm-1
            minWL       = 1e7./[FitLimits{2,1:2}]; % min WL = max cm-1
            ppnmT       = [FitLimits{4,1:2}];

            RefX        = 1e7./RefX;
            MeasY       = flipud(MeasY);
            
            % Remove NaNs from measured data
            MeasX = MeasX(~isnan(MeasY));
            MeasY = MeasY(~isnan(MeasY));
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
%             minRefX  = 345;
%             maxRefX  = 690;

            ppnm = ppnmT;
            wp1 = CWL(1) - 260./ppnm;

            minRefX = minWL(1);
            maxRefX = maxWL(1);

            cutID    = sort(findClosestId2Val(RefX,[minRefX maxRefX]));
            [RefX_cut,srt_id] = sort(RefX(cutID(1):cutID(2)));
            RefY_cut = RefY(cutID(1):cutID(2));
            RefY_cut = RefY_cut(srt_id);

            RefY     = RefY./max(RefY_cut);
            RefY_cut = RefY_cut./max(RefY_cut);
            
        case 3     % UniGE nsTA 
%             minRefX  = 370;
%             maxRefX  = 690;
            ppnm = ppnmT;
            wp1 = CWL(1) - 260./ppnm;

            minRefX = minWL(1);
            maxRefX = maxWL(1);

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
            minRefX  = 830;
            maxRefX  = 1450;

            cutID    = sort(findClosestId2Val(RefX,[minRefX maxRefX]));
            [RefX_cut,srt_id] = sort(RefX(cutID(1):cutID(2)));
            RefY_cut = RefY(cutID(1):cutID(2));
            RefY_cut = RefY_cut(srt_id);

            RefY     = RefY./max(RefY_cut);
            RefY_cut = RefY_cut./max(RefY_cut);
            
            MeasX    = (1:512)';
            MeasY    = MeasY(MeasX);
        case {6,7} % RAL LIFEtime (Absorbance or intensity)
            ppnm     = ppnmT(i); % estimated resolution in pixels/nm
%             dW       = [-570 450]; % min/max relative to CWL to consider for fit
%             dW       = [-400 260]; % min/max relative to CWL to consider for fit
            dW       = [minWL(i) maxWL(i)] - CWL(i);

            wp1      = CWL(i) - 64/ppnm;      % lowest wavelength (nm) = CWL - Npix/2*resolution
            cutID    = sort(findClosestId2Val(RefX,CWL(i)+dW'));
            RefX_cut = RefX(cutID(1):cutID(2));
            RefY_cut = RefY(cutID(1):cutID(2));
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
            % Parameters: 1=centralWL 2=nm/pix 3=Yscale 4=Yshift 5:6=lin+quadratic baseline            
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
%             p0 = [310 1.05  0.6 -0.1 1 1];
            p0 = [wp1       ppnm  0.6 0.1 0.1 0.1];
            LB = [wp1-100   0.5   0.1 -10 -10 -10];
            UB = [wp1+100   2     10   10  10 10 ];
            ftol= 5e-7;
            stol= 5e-7;
        case 3 % UniGE nsTA
%             p0 = [330 1.33  0.58 -0.1 1     1];
            p0 = [wp1       1.33  0.58 -0.1 1     1];
            LB = [wp1-100   0.1   0.1  -1 -Inf -Inf];
            UB = [wp1+100   2     10   1  Inf  Inf ];
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
            ftol= 5e-7;
            stol= 5e-7;
        case {6,7} % RAL LIFEtime
            wp1_L  = CWL(i) - 64/ppnm + dW(1);
            wp1_U  = CWL(i) - 64/ppnm + dW(2);

            p0 = [wp1    ppnm  1.5  -0.1  eps  eps  ];
            LB = [wp1_L  1e-4  0.1  -1    -Inf -Inf ];
            UB = [wp1_U  0.4   10   1     Inf  Inf  ];
            ftol= 5e-7;
            stol= 5e-7;
    end
    
    if doFit == 1
        mb = msgbox('Fit in progress!','Fitting...','help');
        
        options     = optimoptions(@lsqnonlin,...
                    'FunctionTolerance',ftol,...
                    'MaxIterations',5e3,...
                    'MaxFunctionEvaluations',5e4,...
                    'steptolerance',stol,...
                    'display','off',...
                    'typicalX',p0);%,...
%                     'plotfcn','optimplotresnorm_log');
    
        Pfit(i,:) = lsqnonlin(fit_fun,p0,LB,UB,options);
        delete(mb);
    else
        Pfit(i,:)=p0;
    end

    switch CalType
        case {1,4,6,7} % IR Setups UoS, UZH Lab 2, RAL LIFEtime
            MeasX       = 1:Npix(i);
            lam(:,i)    = flipud(MeasX(:)./Pfit(i,2) + Pfit(i,1));
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
        case {1,4,6,7} % IR setups UoS and UZH Lab 2
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