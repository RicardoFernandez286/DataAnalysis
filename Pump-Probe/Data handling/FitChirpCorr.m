function chirpSt = FitChirpCorr(dataStruct,rootfolder,k)
% This routine will fit a dispersion curve by considering a Gaussian IRF
% with a wavelength-dependent t0. These t0's will be fit to Cauchy's
% dispersion equation (generalised to an arbitrary number of terms) to
% give a general time-dependent chirp, which can then be saved to a file.
%
% The 'Automatic' branch uses variable projection (VARPRO): the linear
% amplitudes (Gaussian, derivatives, offset, exponential) are eliminated
% in closed form at each iteration, leaving only (t0, FWHM, tau_exp) as
% nonlinear parameters. trust-region-reflective is used so that bounds
% are genuinely honoured.
%
% A two-pass option is provided. Pass 1 fits each pixel independently
% (cold start). If do2ndPass is true, an intermediate dispersion curve is
% fit to t0(lambda) and a moving-median smoothing of FWHM and tau_exp is
% computed; these provide the initial guess for pass 2, which re-fits
% each pixel and overwrites Pfit.
%
% Pixels with non-finite or all-zero data (masked pump scatter, dead
% detector pixels, etc.) are skipped silently and recorded as NaN in
% Pfit. Both the intermediate and the final dispersion fits filter out
% these pixels.
%
% Ricardo Fernández-Terán / v3.2a / 03.05.2026

warning('off','optimlib:levenbergMarquardt:InfeasibleX0');
warning('off','MATLAB:rankDeficientMatrix');

debugChirp = 0;
%% Read from dataStruct
delays      = dataStruct.delays;
probeAxis   = dataStruct.cmprobe{k};

Ndelays = length(delays);
Npixels = length(probeAxis);

switch dataStruct.rawcorr
    case 'CORRECTED'
        Z = dataStruct.corrdata{k};
    case 'RAW'
        Z = dataStruct.rawsignal{k};
end

chirpSt = [];
%% Hardcoded settings
ChirpEquation     = 'Cauchy';        % 'Cauchy' or 'Shaper'
nCauchyTerms      = 6;               % # of Cauchy terms (>=1): a + b/L^2 + c/L^4 + ...
CauchyLambdaRef   = [];              % Reference wavelength (nm) for parameter scaling. [] uses mean(fitWL).
c0                = 2.99792458e+2;   % Speed of light in nm/fs
plotpercent       = 50;

do2ndPass         = true;            % Run a second VARPRO pass with smoothed inits
smoothWindow      = 10;               % Moving-median window (pixels) for FWHM/tau_exp smoothing

%% Read Probe fit ranges
% probeRanges = [360 380 415 645];
% probeRanges = [0 1000];

%% Select data
% ProbeIdx    = findClosestId2Val(probeAxis,probeRanges);
% fitPixels   = [];
%
% for j=1:round(length(probeRanges)/2)
%     fitPixels = [fitPixels ProbeIdx(2*(j-1)+1):1:ProbeIdx(2*j)];
% end
fitPixels   = 1:Npixels;
fitWL       = probeAxis(fitPixels);
NfitPix     = length(fitPixels);

mode        = questdlg('Select Chirp Correction Method','Chirp Correction Method','Automatic','Manual','Step Function (Auto)','Automatic');

if isempty(mode)
    return
end

switch mode
    case 'Manual'
        % Make a plot of the data
        maxPl = plotpercent./100.*max(abs(Z(:)));
        minPl = -maxPl;
        ctrs  = linspace(minPl,maxPl,40);

        fh = figure(1);
        clf(fh);
        fh.Color = 'w';
        ax = axes('parent',fh);

        % Automatically detect discontinuities in the probe axis (e.g. masked pump scatter)
        % and plot them accordingly
        dProbe          = diff(probeAxis) - mean(diff(probeAxis));
        jump_ID         = dProbe >= 5*mean(diff(probeAxis));
        Zplot           = Z;
        Zplot(:,jump_ID)= NaN;

        contourf(ax,probeAxis,delays,Zplot,ctrs,'EdgeColor','flat');

        colormap(darkb2r(minPl,maxPl,40,2));
        caxis(ax,[minPl,maxPl]);
        colorbar(ax);

        ylim([-20,20]); % Set this preliminary zoom level in case the data contains longer delays

        % Select a region to zoom in
        title(ax,'Define min/max delay region to zoom in','FontWeight','bold');
        ds = [];
        [ds,aborted] = SelectTracesCH(ds,0);
        if aborted ~= 0
            return
        end

        ylim(ax,[min(ds.SelTraces(:,2)) max(ds.SelTraces(:,2))]);

        % Now select points for chirp corr
        title(ax,'Select Data Points for Chirp Correction','FontWeight','bold');
        ylabel(ax,'Delay (ps)','FontWeight','bold');
        xlabel(ax,'Wavelength (nm)','FontWeight','bold');
        ax.FontSize = 16;
        ds = [];
        [ds,aborted] = SelectTracesCH(ds,0);
        if aborted ~= 0
            return
        end

        fitWL       = ds.SelTraces(:,1);
        Pfit(:,1)   = ds.SelTraces(:,2);
        delete(fh);

    case 'Automatic'
        %%% VARPRO: fit Gaussian + (optional) derivatives + (optional) Heaviside-conv-Gaussian
        %%% exponential, eliminating linear amplitudes in closed form.
        %
        % Pfit column layout (kept compatible with the original 8-col scheme):
        %   1 = t0  (ps)
        %   2 = FWHM (ps)
        %   3 = Gaussian amplitude
        %   4 = 1st-derivative amplitude (0 if disabled)
        %   5 = 2nd-derivative amplitude (0 if disabled)
        %   6 = offset
        %   7 = tau_exp (ps)
        %   8 = exponential amplitude

        inclDeriv = menu('Include 1st and 2nd derivatives of the Gaussian IRF?', ...
              'No', '1st Derivative', '2nd Derivative', '1st+2nd Derivatives');
        
        if inclDeriv == 0
            return
        end
        
        switch inclDeriv
            case 1  % 'No'
                derivMask = [false false];
            case 2  % '1st Derivative'
                derivMask = [true  false];
            case 3  % '2nd Derivative'
                derivMask = [false true ];
            case 4  % '1st+2nd Derivatives'
                derivMask = [true  true ];
        end
        useExp = true;

        Pfit = nan(NfitPix, 8);                   % NaN means "pixel skipped"
        Dfit = zeros(Ndelays, NfitPix);
        nf   = ones(NfitPix, 1);

        opts_nl = optimoptions(@lsqnonlin, ...
            'Algorithm','trust-region-reflective', ...
            'Display','off', ...
            'FunctionTolerance',5e-10, ...
            'StepTolerance',5e-10, ...
            'OptimalityTolerance',5e-10, ...
            'MaxIterations',5e2, ...
            'MaxFunctionEvaluations',5e3);

        % ============ PASS 1: cold start at each pixel ============
        wb = waitbar(0,'Pass 1: fitting coherent artefact...');
        figure(wb);

        for j = 1:NfitPix
            i        = fitPixels(j);
            yData    = Z(:,i);

            [~,idM]  = max(abs(yData));   % NaNs ignored; idM=1 if all-NaN/zero

            % Bounds on nonlinear params: [t0, FWHM, tau_exp]
            LB_nl = [delays(idM)-1, 1e-3, 1e-3];
            UB_nl = [delays(idM)+1, 1,    5  ];

            % Cold-start initial guess
            pnl0  = [delays(idM), 0.1, 1];

            [Pfit(j,:), Dfit(:,j), nf(j)] = fitOnePixel(yData, delays, ...
                pnl0, LB_nl, UB_nl, derivMask, useExp, opts_nl);

            waitbar(j/NfitPix,wb,['Pass 1: fitting... (' num2str(j) ' of ' num2str(NfitPix) ')'])

            if debugChirp == 1
                if j == 1
                    fh = figure(2);
                    clf(fh)
                    ax = axes('parent',fh);
                    xlabel('Delay (ps)', 'FontSize',16, 'FontWeight','bold');
                    ylabel('{\Delta}A (mOD)', 'FontSize',16, 'FontWeight','bold')
                end
                plot(ax,delays,Z(:,i),'o')
                hold(ax,'on')
                plot(ax,delays,Dfit(:,j),'-')
                hold(ax,'off')
                title(ax,['Pass 1 - \lambda = ' num2str(probeAxis(i),4) ' nm'])
                xlim(ax,[-1,3]);
                drawnow;
            end
        end
        delete(wb);

        % ============ PASS 2: smoothed initial guesses ============
        if do2ndPass
            % Compute smoothed initial guesses:
            %   t0   from a fitted dispersion curve (Cauchy) or moving median (otherwise)
            %   FWHM from moving-median smoothing
            %   tau  from moving-median smoothing
            [t0_init, fwhm_init, tau_init] = computeSecondPassInits(Pfit, fitWL, ...
                ChirpEquation, nCauchyTerms, CauchyLambdaRef, smoothWindow);

            wb = waitbar(0,'Pass 2: refining with smoothed initial guesses...');

            for j = 1:NfitPix
                i      = fitPixels(j);
                yData  = Z(:,i);

                [~,idM]  = max(abs(yData));

                % Fall back to peak-of-data if the dispersion init is non-finite
                if isfinite(t0_init(j))
                    t0_guess = t0_init(j);
                else
                    t0_guess = delays(idM);
                end

                % Bounds anchored on the t0 estimate
                LB_nl = [t0_guess-1, 1e-3, 1e-3];
                UB_nl = [t0_guess+1, 1,    5  ];

                % Initial guess from the smoothed neighbourhood
                %   max(NaN, 1e-3) = 1e-3 in MATLAB, so this is NaN-safe.
                pnl0  = [t0_guess, max(fwhm_init(j),1e-3), max(tau_init(j),1e-3)];

                [Pfit(j,:), Dfit(:,j), nf(j)] = fitOnePixel(yData, delays, ...
                    pnl0, LB_nl, UB_nl, derivMask, useExp, opts_nl);

                waitbar(j/NfitPix,wb,['Pass 2: refining... (' num2str(j) ' of ' num2str(NfitPix) ')'])

                if debugChirp == 1
                    if j == 1
                        fh = figure(3);
                        clf(fh)
                        ax = axes('parent',fh);
                        xlabel('Delay (ps)', 'FontSize',16, 'FontWeight','bold');
                        ylabel('{\Delta}A (mOD)', 'FontSize',16, 'FontWeight','bold')
                    end
                    plot(ax,delays,Z(:,i),'o')
                    hold(ax,'on')
                    plot(ax,delays,Dfit(:,j),'-')
                    hold(ax,'off')
                    title(ax,['Pass 2 - \lambda = ' num2str(probeAxis(i),4) ' nm'])
                    xlim(ax,[-1,3]);
                    drawnow;
                end
            end
            delete(wb);
        end

        % Report skipped pixels (data NaN/zero or fit failed)
        nSkipped = sum(~isfinite(Pfit(:,1)));
        if nSkipped > 0
            fprintf('FitChirpCorr: %d of %d pixels skipped (non-finite/empty data or fit failure).\n', ...
                nSkipped, NfitPix);
        end

    case 'Step Function (Auto)'
        %%% Fit a Gaussian response + step and take the WL-dependent t0 parameter to do the chirp fit
        % 1D Gaussian with 1st and 2nd derivative
        % Parameters: 1=t0  2=FWHM  3=Amplitude  4=Amp 1st deriv  5=Amp 2nd deriv   6 = Offset
        G0 = @(p,t) exp(-log(2).*((t-p(1))./p(2)).^2);
        G1 = @(p,t) -2.*log(2)./(p(2).^2).*(t-p(1)).*G0(p,t);
        G2 = @(p,t) 2.*log(2).*(-p(2).^2 + 2.*log(2).*(t-p(1)).^2)./(p(2).^4).*G0(p,t);
        HS = @(p,t) p(7)*(1+erf((t-p(1))./(sqrt(2).*p(2)/(2*sqrt(2*log(2))))));

        ArtFitStep = @(p,t) p(3).*G0(p,t) + p(4).*G1(p,t) + p(5).*G2(p,t) + p(6) + HS(p,t);

        P0  = [0.1  0.1   1     1    1    1     1  ];
        LB  = [-5   0     -50   -50  -50  -50   -50];
        UB  = [5    1     50    50   50   50    50 ];

        % Trim the data to consider only early delays
        %   Typically [-1,3] ps works well for TA/TRIR
        tmin = -1;
        tmax = 3;
        fitdelID = delays >= tmin & delays <= tmax;
        fit_delays = delays(fitdelID);
        Ndelays = length(fit_delays);
        Zfit = Z(fitdelID,:);

        % Prepare an array where we will store the results
        Pfit = zeros(NfitPix,length(P0));
        Dfit = zeros(Ndelays,NfitPix);
        % Define fit options
        options     = optimoptions(@lsqcurvefit,...
                        'FunctionTolerance',5e-10,...
                        'MaxIterations',1e4,...
                        'MaxFunctionEvaluations',1e4,...
                        'steptolerance',5e-10,...
                        'OptimalityTolerance',5e-10,...
                        'display','off');

        inclDeriv = menu('Include 1st and 2nd derivatives of the Gaussian IRF?', ...
              'No', '1st Derivative', '2nd Derivative', '1st+2nd Derivatives');

        if inclDeriv == 0
            figure(app.DataAnalysisGUI_UIFigure);
            return
        end

        % Do the fits
        wb = waitbar(0);
        for j=1:NfitPix
            i=fitPixels(j);
            [P0(3),idM] = max(abs(diff(abs(Zfit(:,i)))));
            P0(1)   = fit_delays(idM);
            UB(1)   = P0(1) + 1;
            LB(1)   = P0(1) - 1;
            UB(3:7) = 10*P0(3);
            LB(3:7) = -10*P0(3);

            % Do NOT include derivative terms (?)
            switch inclDeriv
                case 1  % 'No'
                    LB(4:5) = 0;
                    UB(4:5) = 0;
                case 2  % '1st Derivative'
                    UB(5)   = 0;
                    LB(5)   = 0;
                case 3  % '2nd Derivative'
                    UB(4)   = 0;
                    LB(4)   = 0;
            end
            Pfit(j,:) = lsqcurvefit(ArtFitStep,P0,fit_delays,Zfit(:,i),LB,UB,options);
            Dfit(:,j) = ArtFitStep(Pfit(j,:),fit_delays);
            waitbar(j/NfitPix,wb,['Fitting chirp correction... (' num2str(j) ' of ' num2str(NfitPix) ')'])
        end
        delete(wb);

end

%% Diagnostic plot for fit
% j=80;
% figure(5)
% plot(delays,Dfit(:,j),'-','LineWidth',2)
% hold on
% plot(delays,Z(:,fitPixels(j)),'o');
% hold off;
%
% Pfit(j,:)

%% Fit Chirp vs Lambda

switch ChirpEquation
    case 'Cauchy'
        % Generalised Cauchy dispersion truncated at nCauchyTerms:
        %   t0(L) = p(1) + sum_{n=1}^{nCauchyTerms-1} ( lambda_ref^(2n) * p(n+1) ) / L^(2n)
        %
        % The lambda_ref^(2n) prefactors render all higher-order
        % coefficients dimensionless and of comparable magnitude, which
        % keeps the least-squares problem well-conditioned regardless of
        % nCauchyTerms.
        nTerms      = max(1, round(nCauchyTerms));
        if isempty(CauchyLambdaRef)
            lambda_ref = mean(fitWL,'omitnan');
        else
            lambda_ref = CauchyLambdaRef;
        end
        expVec      = 2*(0:nTerms-1);                 % [0 2 4 6 ...]
        scaleVec    = lambda_ref.^expVec;             % [1 lambda_ref^2 lambda_ref^4 ...]

        chirpFun    = @(p,L) reshape( ...
                        sum( (scaleVec .* reshape(p,1,[])) ./ (L(:).^expVec), 2 ), ...
                        size(L) );

        % Initial guess and bounds (parameters are O(1) thanks to the scaling)
        C0 = [3,    zeros(1, nTerms-1)];
        UB = [5,    0*1e2*ones(1, nTerms-1)];
        LB = [-5,  -1e2*ones(1, nTerms-1)];
    case 'Shaper'
        % Fit dispersion using a Taylor expansion centred at w0 --- ONLY FOR PULSE SHAPER CALIBRATION
        % Parameters: 1=a 2=b 3=c; (4=d)  L = Wavelength in nm   W0 = central wavelength

        % Need to work in THz and fs
        Pfit(:,1) = Pfit(:,1).*1000;
        delays    = delays.*1000;

        % % If probe = cm-1
        % probeAxis = 1e7./probeAxis;
        % fitWL     = 1e7./fitWL;

        CWL       = mean([min(fitWL) max(fitWL)]);

        probeAxis = 1e3*c0./probeAxis;
        fitWL     = 1e3*c0./fitWL;

        % Expand around the central wavelength
        nu0     = c0./CWL;
        SpPhase = @(p,nu) p(2).*(2.*pi).*(nu-nu0) + 1/2*p(3).*(2.*pi).^2.*(nu-nu0).^2 + 1/6*p(4).*(2.*pi).^3.*(nu-nu0).^3 + 1/24*p(5).*(2.*pi).^4.*(nu-nu0).^4;
        chirpFun= @(p,nu) p(1) + 1./nu.*SpPhase(p,nu);

%         C0 = [1.3e3 125   -25    20    -20 ];
        C0 = [1.3e3  eps   -1.5   -1.5   0  ];
        UB = [1e6   1E10    1e10    0   0 ];
        LB = [-1e6  -1E10  -1e10   -0  0];
end

%         options     = optimoptions(@lsqcurvefit,...
%                         'FunctionTolerance',1e-10,...
%                         'MaxIterations',5e4,...
%                         'MaxFunctionEvaluations',5e4,...
%                         'steptolerance',1e-10,...
%                         'display','off');
%                         'typicalX',C0,...

        opt2            = optimset(@lsqcurvefit);
        opt2.TolFun     = 1e-10;
        opt2.MaxIter    = 5e4;
        opt2.MaxFunEvals= 5e4;
        opt2.TolX       = 1e-10;
        opt2.Display    = 'off';

% Filter out non-finite t0 values (from skipped pixels) before the dispersion fit
valid_disp = isfinite(Pfit(:,1));
% Cfit    = lsqcurvefit(chirpFun,C0,fitWL,Pfit(:,1),LB,UB,options);
Cfit    = robustlsqcurvefit(chirpFun,C0,fitWL(valid_disp),Pfit(valid_disp,1),LB,UB,'bisquare',opt2);
% Cfit = C0;

% Plot Everything
%%% Plot fitted chirp data
WHAT    = Z;
maxPl   = max(abs(WHAT(:)));
minPl   = -maxPl;
ctrs    = linspace(minPl,maxPl,40);

fh      = figure(1);
clf(fh);
fh.Color= 'w';

ax      = axes('parent',fh);

% Automatically detect discontinuities in the probe axis (e.g. masked pump scatter)
% and plot them accordingly
dProbe          = diff(probeAxis) - mean(diff(probeAxis));
jump_ID         = dProbe >= 5*mean(diff(probeAxis));
WHAT(:,jump_ID) = NaN;

contourf(ax,probeAxis,delays,WHAT,ctrs,'EdgeColor','flat');
colormap(darkb2r(minPl,maxPl,40,2));
caxis([minPl,maxPl]);
colorbar;

hold(ax,'on');
plot(ax,fitWL,Pfit(:,1),'o','Color',0.8*[0 1 1])
plot(ax,probeAxis,chirpFun(Cfit,probeAxis),'y','LineWidth',3);
hold(ax,'off');

switch ChirpEquation
    case 'Cauchy'
        ylim(ax,[mean(Pfit(:,1),'omitnan')-2 mean(Pfit(:,1),'omitnan')+2]) % time is in ps
        title(ax,'Fitted Dispersion Curve','FontWeight','bold');
        ylabel(ax,'Delay (ps)','FontWeight','bold');
        xlabel(ax,'Wavelength (nm)','FontWeight','bold');
    case 'Shaper'
        ylim(ax,[min(Pfit(:,1))-250 max(Pfit(:,1))+250]) % time is in fs
        title(ax,'Fitted Dispersion Curve','FontWeight','bold');
        ylabel(ax,'Delay (fs)','FontWeight','bold');
        xlabel(ax,'Frequency (THz)','FontWeight','bold');
        ax.XDir='reverse';
end
ax.FontSize = 16;

%%% Plot dispersion fit
fh = figure(2);
clf(fh);
fh.Color = 'w';
ax = axes('parent',fh);

plot(ax,fitWL,Pfit(:,1),'xk')
hold(ax,'on')
plot(ax,probeAxis,chirpFun(Cfit,probeAxis),'-r','LineWidth',1.5)
hold(ax,'off')
axis(ax,'tight');
title(ax,'Fitted Dispersion Curve','FontWeight','bold');
ylabel(ax,'Delay (ps)','FontWeight','bold');
axis(ax,'auto y');
box(ax,'on');

fh.Position(3:4) = [640 500];
fh.Position(2) = 100;


ax2 = axes('parent',fh);

ax.Position = [0.15 0.46 0.7750 0.4612];
ax2.Position= [0.15 0.15 0.7750 0.2528];

% Multiply residuals (delta t_0) to adjust time scales (always plot in fs)
switch ChirpEquation
    case 'Cauchy'
        factorRes   = 1000; % ps to fs
    case 'Shaper'
        factorRes   = 1; % fs
end

res = factorRes.*(Pfit(:,1)-chirpFun(Cfit,fitWL));
plot(ax2,fitWL,res,'xk')
yline(ax2,0);
ax.XTickLabel=[];

switch ChirpEquation
    case 'Cauchy'
        ylim(ax,[min(Pfit(:,1))-0.25 max(Pfit(:,1))+0.25])
        xlabel(ax2,'Wavelength (nm)','FontWeight','bold');
        ylabel(ax2,'{\Delta}{t_{0}}^{res} (fs)','FontWeight','bold');
        ylim(ax2,[-50,50]); % Always plot +/- 50 fs delta t_0
        xlim(ax,[min(probeAxis) max(probeAxis)])
    case 'Shaper'
%         ylim(ax,[min(Pfit(:,1))-250 max(Pfit(:,1))+250])
        title(ax,'Fitted Dispersion Curve','FontWeight','bold');
        ylabel(ax,'Delay (fs)','FontWeight','bold');
        xlabel(ax2,'Frequency (THz)','FontWeight','bold');
        xline(ax,nu0);
        ax2.XDir='reverse';
        ax.XDir='reverse';
%         ylim(ax,[-pi,pi]);
%         yline(ax,Cfit(1));
end

title(ax,'Fitted Dispersion Curve','FontWeight','bold');

box(ax2,'on');
ax.FontSize = 16;
ax2.FontSize = 16;
ax.TickLength = [0.02 0.02];
ax2.TickLength = [0.02 0.02];
linkaxes([ax,ax2],'x');

switch mode
    case {'Automatic','Step Function (Auto)'}
        %%% Plot wavelength-dependent FWHM
        % Pfit(:,2) is FWHM (in ps) for the 'Automatic' branch but sigma
        % (in ps) for the 'Step Function (Auto)' branch — apply the
        % conversion factor so that the displayed quantity is always FWHM.
        switch mode
            case 'Automatic'
                fwhmFactor = 1;                     % already FWHM
            case 'Step Function (Auto)'
                fwhmFactor = 2*sqrt(2*log(2));      % sigma -> FWHM
        end

        fh = figure(3);
        clf(fh);
        fh.Color = 'w';
        ax = axes('parent',fh);

        plot(ax,fitWL,Pfit(:,2)*fwhmFactor*1000,'xk')
        axis(ax,'tight');
        xlabel(ax,'Wavelength (nm)','FontWeight','bold');
        ylabel(ax,'IRF FWHM (fs)','FontWeight','bold');
        ylim(ax,[0 500]);
        xlim(ax,[min(probeAxis) max(probeAxis)]);
        box(ax,'on');

        meanIRF = mean(Pfit(:,2),'omitnan')*fwhmFactor*1000;

        yline(ax,meanIRF,'--','LineWidth',2,'Color','b');

        title(ax,['Fitted IRF FWHM - Avg: ' num2str(meanIRF,'%.1f') ' fs'],'FontWeight','bold');

        fh.Position(3:4) = [640 370];
        ax.FontSize = 16;
        ax.TickLength = [0.02 0.02];
end

f = msgbox('Fit done! Please verify the results, then click OK.');
uiwait(f);

keepFit = questdlg('Are you happy with these results?','Happy?','Yes, keep them!','Nope :(','Nope :(');
switch keepFit
    case 'Nope :('
        return
    case 'Yes, keep them!'
        % Save the chirp correction .mat file
        path_save = uigetdir(rootfolder,'Save Chirp Correction Parameters to...');
        fn_save = ['chirpCorr_' datestr(datetime('now'),'yyyymmdd-HHMMSS') '.mat'];
        try
            save([path_save filesep fn_save],'chirpFun','Cfit','fitPixels','fitWL','Pfit');
            helpdlg(["Chirp correction file successfully saved to: "; [path_save filesep fn_save]],'Chirp Correction Saved!');
        catch ME
            errordlg('Error saving chirp correction file','Error saving file');
        end

        % Copy variables to output structure
        chirpSt.chirpFun    = chirpFun;
        chirpSt.Cfit        = Cfit;
        chirpSt.fitPixels   = fitPixels;
        chirpSt.fitWL       = fitWL;
        chirpSt.Pfit        = Pfit;
        switch mode
            case 'Automatic'
                chirpSt.meanIRF     = meanIRF;
            otherwise
                chirpSt.meanIRF     = NaN;
        end
end

end

%% ===== Local helpers =====================================================

function [pfit_row, dfit_col, nf_val] = fitOnePixel(yData, delays, ...
    pnl0, LB_nl, UB_nl, derivMask, useExp, opts_nl)
% Run a single VARPRO pixel fit and pack the result into the 8-column
% Pfit row layout.
%
%   pfit_row : 1x8 row vector (NaN if pixel skipped or fit failed)
%   dfit_col : Ndelays x 1 reconstructed model (zeros if skipped)
%   nf_val   : normalisation factor (max(|yData|)) used during the fit
%
% Pixels whose data are all-NaN, all-zero, or contain Inf are skipped
% silently. Fit failures are also caught so that one bad pixel does not
% abort the whole loop.

    pfit_row = nan(1, 8);
    dfit_col = zeros(numel(delays), 1);
    nf_val   = 1;

    % --- Data sanity checks ---
    if ~all(isfinite(yData))
        return;   % NaN/Inf in column (e.g. masked pump scatter, dead detector)
    end
    nf_test = max(abs(yData));
    if ~isfinite(nf_test) || nf_test == 0
        return;   % Constant-zero column
    end

    nf_val  = nf_test;
    yData_n = yData ./ nf_val;

    % Make sure the initial guess is feasible
    pnl0 = min(max(pnl0, LB_nl), UB_nl);
    if any(~isfinite(pnl0))
        return;
    end

    % --- Run VARPRO ---
    try
        pnl = lsqnonlin(@(p) varpro_resid(p, delays, yData_n, derivMask, useExp), ...
                        pnl0, LB_nl, UB_nl, opts_nl);
    catch
        % Optimiser failed: keep pixel as NaN
        return;
    end

    [~, A, amps] = varpro_resid(pnl, delays, yData_n, derivMask, useExp);

    pfit_row     = zeros(1, 8);
    pfit_row(1)  = pnl(1);    % t0
    pfit_row(2)  = pnl(2);    % FWHM
    pfit_row(7)  = pnl(3);    % tau_exp

    ampIdx       = 1;
    pfit_row(3)  = amps(ampIdx); ampIdx = ampIdx + 1;            % G0 amp
    if derivMask(1)
        pfit_row(4) = amps(ampIdx); ampIdx = ampIdx + 1;         % G1 amp
    end
    if derivMask(2)
        pfit_row(5) = amps(ampIdx); ampIdx = ampIdx + 1;         % G2 amp
    end
    pfit_row(6)  = amps(ampIdx); ampIdx = ampIdx + 1;            % offset
    if useExp
        pfit_row(8) = amps(ampIdx);                              % exp amp
    end

    dfit_col = nf_val .* (A * amps);
end

function [t0_init, fwhm_init, tau_init] = computeSecondPassInits(Pfit, fitWL, ...
    ChirpEquation, nCauchyTerms, CauchyLambdaRef, smoothWindow)
% Build smoothed initial guesses for the 2nd VARPRO pass.
%   t0_init    : from a fitted dispersion curve (Cauchy) or movmedian (else)
%   fwhm_init  : moving-median smoothed FWHM
%   tau_init   : moving-median smoothed tau_exp
%
% NaN entries in Pfit are excluded from medians (via 'omitnan') and from
% the intermediate dispersion fit. Output vectors share the shape of fitWL.

    fwhm_init = movmedian(Pfit(:,2), smoothWindow, 'omitnan');
    tau_init  = movmedian(Pfit(:,7), smoothWindow, 'omitnan');

    valid = isfinite(Pfit(:,1));

    switch ChirpEquation
        case 'Cauchy'
            nTerms = max(1, round(nCauchyTerms));
            if isempty(CauchyLambdaRef)
                lambda_ref = mean(fitWL,'omitnan');
            else
                lambda_ref = CauchyLambdaRef;
            end
            expVec   = 2*(0:nTerms-1);
            scaleVec = lambda_ref.^expVec;

            chirpFun_local = @(p,L) reshape( ...
                sum( (scaleVec .* reshape(p,1,[])) ./ (L(:).^expVec), 2 ), ...
                size(L) );

            C0 = [3,   zeros(1, nTerms-1)];
            UB = [5,   1e2*ones(1, nTerms-1)];
            LB = [-5, -1e2*ones(1, nTerms-1)];

            opt2             = optimset(@lsqcurvefit);
            opt2.TolFun      = 1e-10;
            opt2.MaxIter     = 5e4;
            opt2.MaxFunEvals = 5e4;
            opt2.TolX        = 1e-10;
            opt2.Display     = 'off';

            if nnz(valid) >= nTerms
                try
                    Cfit_loc = robustlsqcurvefit(chirpFun_local, C0, ...
                        fitWL(valid), Pfit(valid,1), LB, UB, 'bisquare', opt2);
                    t0_init = chirpFun_local(Cfit_loc, fitWL);
                catch
                    % Robust fit failed; fall back to median smoothing
                    t0_init = movmedian(Pfit(:,1), smoothWindow, 'omitnan');
                end
            else
                % Not enough valid points for a parametric fit
                t0_init = movmedian(Pfit(:,1), smoothWindow, 'omitnan');
            end

        otherwise
            % Fallback: median-smooth t0 directly
            t0_init = movmedian(Pfit(:,1), smoothWindow, 'omitnan');
    end
end

function [resid, A, amps] = varpro_resid(pnl, t, y, derivMask, useExp)
% VARPRO residual for the coherent-artefact model.
%   pnl       = [t0, FWHM, tau_exp]
%   derivMask = logical 1x2 -> [include G1, include G2]
%   useExp    = logical -> include erf-broadened exponential term

    t  = t(:);
    y  = y(:);
    t0 = pnl(1);
    sg = pnl(2)/(2*sqrt(2*log(2)));   % FWHM -> sigma
    tau= pnl(3);

    G0 = exp(-log(2).*((t-t0)./sg).^2);
    A  = G0;

    if derivMask(1)
        G1 = -2.*log(2)./(sg.^2).*(t-t0).*G0;
        A  = [A G1];
    end
    if derivMask(2)
        G2 = 2.*log(2).*(-sg.^2 + 2.*log(2).*(t-t0).^2)./(sg.^4).*G0;
        A  = [A G2];
    end

    A = [A ones(numel(t),1)];   % offset

    if useExp
        % Numerically stable, piecewise-equivalent erf-broadened exponential.
        %   For (t-t0) >= 0 : standard form (exp prefactor decays).
        %   For (t-t0) <  0 : erfcx-based form (Gaussian prefactor decays).
        % Both forms are mathematically identical and continuous at t = t0.
        u    = t - t0;
        sgsq = sg.^2;
        late = u >= 0;
        E1   = zeros(numel(u),1);

        if any(late)
            ul        = u(late);
            argE      = -ul./tau + 0.5.*sgsq./tau.^2;
            argB      = (ul - sgsq./tau)./(sg.*sqrt(2));
            E1(late)  = 0.5.*exp(argE).*(1+erf(argB));
        end
        if any(~late)
            ue           = u(~late);
            prefac       = exp(-ue.^2./(2.*sgsq));
            arg_x        = sg./(tau.*sqrt(2)) - ue./(sg.*sqrt(2));
            E1(~late)    = 0.5.*prefac.*erfcx(arg_x);
        end

        A = [A E1];
    end

    amps  = A \ y;
    resid = y - A*amps;
end

function [dataStruct,aborted] = SelectTracesCH(dataStruct,doSort,varargin)
% Output is a sorted (or not) vector, dataStruct.SelTraces, of the selected (X,Y) points via interactive input through the mouse.
%
% USAGE:
% Left-click to add points (right-click to remove last)
%
% Press Return key when done.
%
% Press middle button of the mouse to exit
%    (will give an error in the caller routine if nothing was selected)
%
% Ricardo Fernández-Terán
% v2.1b - 25.11.2021

% Clear the output variable (in case it already existed)
dataStruct.SelTraces = [];

warning('off','MATLAB:ginput:FigureDeletionPause');
warning('off','MATLAB:ginput:Interrupted');

if nargin < 2
    doSort = 1;
end

i=0;
hold on
while 1 % Repeat the loop until Return is pressed
    [xp,yp,button] = ginput(1);
    if isequal(button,1)
        i=i+1;
        plot(xp,yp,'o','Linewidth',1.5,'Color','g')
        dataStruct.SelTraces(i,:)=[xp yp];
    elseif isequal(button,3) && i>0
        % Remove last point when right clicking
        plot(dataStruct.SelTraces(i,1),dataStruct.SelTraces(i,2),'x','Linewidth',1.5,'Color','r')
        i=i-1;
    elseif isequal(button,2)
        aborted = 1;
        break
    elseif isempty(xp)
        aborted = 0;
        break
    end
end

dataStruct.SelTraces = dataStruct.SelTraces(1:i,:);

if doSort == 1
    dataStruct.SelTraces = sort(dataStruct.SelTraces,1);
end
hold off
end