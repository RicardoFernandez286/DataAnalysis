function [cm,lam] = CalibrateWLAxis(CAL_data)
% Subroutine to calibrate the wavelength axis of a spectrometer, using a centering and scaling approach.
% Will output the fitted cm-1 and nm axes
% The number of columns of cm and lam = Ndet (given by the size of the AbsSolvent cell in the CAL_data structure).
% Ricardo Fernandez-Teran / 09.11.2021 / v1.0a

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
    if CalType == 1 
        Gratings    = CAL_data.Gratings;
        CWL         = CAL_data.CWL;
        RefX        = 1e7./RefX;
        MeasY       = flipud(MeasY);
    end

    %% Cut off relevant part of reference spectrum
    switch CalType
        case 1 % UoS 2DIR/TRIR
            switch Gratings(i)
                case 0 
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
            minRefX  = 340;
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
    end

    %% Define Fit Functions
    ref_fun     = griddedInterpolant(RefX_cut,RefY_cut);
    meas_fun    = griddedInterpolant(MeasX,MeasY);

    % Parameters: 1=offset 2=nm/pix 3=Yscale 4=Yshift            
    plot_fun    = @(p) meas_fun((RefX_cut-p(1)).*p(2)).*p(3)+p(4)+1e-6.*p(5).*(RefX_cut-p(1))+1e-6.*p(6).*(RefX_cut-p(1)).^2;
    fit_fun     = @(p) (plot_fun(p) - ref_fun(RefX_cut)).*1e4;
    %% Calibration Fit

    switch CalType
        case 1 % UoS TRIR/2DIR
            p0 = [wp1  ppnm  1    -0.1  eps eps];
            LB = [1e3  1e-4  0.1  -1    -1  -1 ];
            UB = [1e4  0.25  10   1     1   1  ];
            ftol= 5e-7;
            stol= 5e-7;
        case 2 % UniGE TA
            p0 = [310 1.05  0.6 -0.1 eps eps];
            LB = [100 0.1   0.1 -1   -1 -1];
            UB = [400 2     10  1    1 1];
            ftol= 5e-7;
            stol= 5e-7;
        case 3 % UniGE nsTA
            p0 = [330 1.33  0.58 -0.1 eps eps ];
            LB = [100 0.1   0.1  -1   -1 -1];
            UB = [400 2     10   1    1 1];
            ftol= 1e-6;
            stol= 1e-6;
    end
    
    mb = msgbox('Fit in progress!','Fitting...','help');
    
    options     = optimoptions(@lsqnonlin,...
                'FunctionTolerance',ftol,...
                'MaxIterations',5e3,...
                'MaxFunctionEvaluations',5e4,...
                'steptolerance',stol,...
                'display','off',...
                'typicalX',p0);%,...
%                 'plotfcn','optimplotresnorm_log');
    
    Pfit(i,:) = lsqnonlin(fit_fun,p0,LB,UB,options);
    
    delete(mb);
    % Pfit(i,:)=p0;

    if CalType ~= 1
        lam(:,i)    = MeasX/Pfit(i,2) + Pfit(i,1); 
    else
        lam(:,i)    = flipud(MeasX/Pfit(i,2) + Pfit(i,1));
    end
    cm(:,i)         = 1e7./lam(:,i);

    %% Plot Fit Results
    fh{i} = figure(i);
    fh{i}.Color = 'w';
    fh{i}.Position(3:4) = [600 350];
    clf(fh{i});
    ax = axes('parent',fh{i});

    if i>1
        fh{i}.Position(2) = fh{i-1}.Position(2) - fh{i-1}.Position(4) - 100;
    end

    switch CalType
        case 1
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
        
        case {2,3}
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