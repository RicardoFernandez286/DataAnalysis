function chirpSt = FitChirpCorr(dataStruct,rootfolder,k)
% This routine will fit a dispersion curve by considering a Gaussian IRF
% with a wavelength-dependent t0. These t0's will be fit to Cauchy's
% dispersion equation to give a general time-dependent chirp, which can
% then be saved to a file.
%
% Ricardo Fernández-Terán / v1.0a / 22.11.2021

warning('off','optimlib:levenbergMarquardt:InfeasibleX0');

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
ChirpEquation   = 'Cauchy'; % 'Cauchy' or 'Shaper';
c0              = 2.99792458e+2;    % Speed of light in nm/fs
plotpercent     = 50;
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
        %%% Fit a Gaussian response and take the WL-dependent t0 parameter to do the chirp fit
        % 1D Gaussian with 1st and 2nd derivative
        % Parameters: 1=t0  2=FWHM  3=Amplitude  4=Amp 1st deriv  5=Amp 2nd deriv  6=Offset  7=ExpTau 8=ExpAmp
        %
        G0 = @(p,t) exp(-log(2).*((t-p(1))./p(2)).^2);
        G1 = @(p,t) -2.*log(2)./(p(2).^2).*(t-p(1)).*G0(p,t);
        G2 = @(p,t) 2.*log(2).*(-p(2).^2 + 2.*log(2).*(t-p(1)).^2)./(p(2).^4).*G0(p,t);
        
        s  = @(p) p(2)/(2*sqrt(2*log(2)));
        E1 = @(p,t) p(8)*0.5*exp(-1./p(7).*(t-p(1) - 0.5.*1/p(7).*s(p).^2)).*(1+erf((t-p(1)-1./p(7).*s(p).^2)./(s(p).*sqrt(2))));

        %%% Fit Gaussian + derivatives
        % ArtFit = @(p,t) p(3).*G0(p,t) + p(4).*G1(p,t) + p(5).*G2(p,t) + p(6);
        % P0  = [0.1  0.1   1     0.1   0.1  0.1     ];
        % LB  = [-5   0     -50   -50  -50  -50   ];
        % UB  = [5    1     50    50   50   50    ];
        
        %% Fit Gaussian + derivatives + exponential
        ArtFit = @(p,t) p(3).*G0(p,t) + p(4).*G1(p,t) + p(5).*G2(p,t) + p(6) + E1(p,t);
        P0  = [0.1  0.1   1     0.1   0.1  0.1  2 1];
        LB  = [-5   0     -10   -10  -10  -10   0 -10];
        UB  = [5    1     10    10   10   10    5 10];

        % % Fit ExpConvGauss only
        % Parameters: 1=t0  2=FWHM  3=Offset 4=ExpTau 5=ExpAmp
        % ArtFit= @(p,t) p(3) + p(5)*0.5*exp(-1./p(4).*(t-p(1) - 0.5.*1/p(4).*s(p).^2)).*(1+erf((t-p(1)-1./p(4).*s(p).^2)./(s(p).*sqrt(2))));
        % P0 = [1 0.1 0.1 0.5 1];
        % UB = [5   1   10  5   10];
        % LB = [-5  0   -10 0   -10];
        
        doFit = 1;

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
                        'Algorithm','levenberg-marquardt',...
                        'display','off');
        
        qf = uifigure;
        qf.Position(3:4) = [650 170];
        inclDeriv = uiconfirm(qf,'Include 1st and 2nd derivatives of the Gaussian IRF?','Artifact Fit Options','Options',{'No';'1st Derivative';'2nd Derivative';'1st+2nd Derivatives'},'DefaultOption',1,'Icon','question');
        delete(qf);

        if strcmp(inclDeriv,'Cancel')
            figure(app.DataAnalysisGUI_UIFigure);
            return
        end

        % Do the fits
        wb = waitbar(0);     
        
        for j=1:NfitPix
            i=fitPixels(j);
            
            % Normalise the data and store norm to reconstruct from fit
            yData = Z(:,i);
            nf(j) = max(abs(yData));
            yData_n = yData./nf(j);

            [~,idM] = max(abs(yData_n));
            P0(3)   = yData_n(idM);
            P0(1)   = delays(idM);
            UB(1)   = P0(1) + 1;
            LB(1)   = P0(1) - 1;
            % UB(3:6) = 10*P0(3);
            % LB(3:6) = -10*P0(3); 
            
            % Do/do not include derivative terms
            switch inclDeriv
                case 'No'
                    LB(4:5) = 0;
                    UB(4:5) = 0;
                case '1st Derivative'
                    UB(5)   = 0;
                    LB(5)   = 0;
                case '2nd Derivative'
                    UB(4)   = 0;
                    LB(4)   = 0;
            end


            
            if doFit == 1
                Pfit(j,:) = lsqcurvefit(ArtFit,P0,delays,yData_n,LB,UB,options);
            else
                Pfit(j,:) = P0;
            end

            Dfit(:,j) = nf(j).*ArtFit(Pfit(j,:),delays);
            waitbar(j/NfitPix,wb,['Fitting chirp correction... (' num2str(j) ' of ' num2str(NfitPix) ')'])
            
            % % Plot the single pixel fits
            % fh=figure(1);
            % clf(fh)
            % ax = axes('parent',fh);
            % hold(ax,'on')
            % plot(delays,Z(:,i),'o')
            % plot(delays,Dfit(:,j),'-')
            % hold(ax,'off')
            % drawnow;
        end
        delete(wb);
        %%
        % return
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
        
        qf = uifigure;
        qf.Position(3:4) = [650 170];
        inclDeriv = uiconfirm(qf,'Include 1st and 2nd derivatives of the Gaussian IRF?','Artifact Fit Options','Options',{'No';'1st Derivative';'2nd Derivative';'1st+2nd Derivatives'},'DefaultOption',1,'Icon','question');
        delete(qf);

        if strcmp(inclDeriv,'Cancel')
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
                case 'No'
                    LB(4:5) = 0;
                    UB(4:5) = 0;
                case '1st Derivative'
                    UB(5)   = 0;
                    LB(5)   = 0;
                case '2nd Derivative'
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
        % Fit dispersion using Cauchy's Equation (truncated on the 3rd/4th term)
        % Parameters: 1=a 2=b 3=c; (4=d)  L = Wavelength in nm
        chirpFun = @(p,L) p(1) + 1e8.*p(2)./(L.^2) + 1e12.*p(3)./(L.^4);% + 1e17.*p(4)./(L.^6);

        C0 = [3     1       1    ];%    1    ];
        UB = [5     1e3     1e3  ];%    1e7  ];
        LB = [-5    -1e3    -1e3 ];%    -1e7 ];
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

% Cfit    = lsqcurvefit(chirpFun,C0,fitWL,Pfit(:,1),LB,UB,options);
% Cfit    = robustlsqcurvefit(chirpFun,C0,fitWL,Pfit(:,1),LB,UB,'bisquare',opt2);
Cfit = C0;

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
        ylim(ax,[mean(Pfit(:,1))-2 mean(Pfit(:,1))+2]) % time is in ps
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
        fh = figure(3);
        clf(fh);
        fh.Color = 'w';
        ax = axes('parent',fh);

        plot(ax,fitWL,Pfit(:,2)*1000,'xk')
        axis(ax,'tight');
        xlabel(ax,'Wavelength (nm)','FontWeight','bold');
        ylabel(ax,'IRF FWHM (fs)','FontWeight','bold');
        ylim(ax,[0 250]);
        xlim(ax,[min(probeAxis) max(probeAxis)]);
        box(ax,'on');

        meanIRF = mean(Pfit(:,2))*1000;

        yline(ax,meanIRF,'--','LineWidth',2,'Color','b');

        title(ax,['Fitted IRF FWHM - Avg: ' num2str(meanIRF,'%.1f') ' fs'],'FontWeight','bold');

        fh.Position(3:4) = [640 370];
        ax.FontSize = 16;
        ax.TickLength = [0.02 0.02];
end

return

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