function chirpSt = FitChirpCorr(dataStruct,rootfolder)
% This routine will fit a dispersion curve by considering a Gaussian IRF
% with a wavelength-dependent t0. These t0's will be fit to Cauchy's
% dispersion equation to give a general time-dependent chirp, which can
% then be saved to a file.
%
% Ricardo Fernández-Terán / v1.0a / 22.11.2021

%% Read from dataStruct
delays      = dataStruct.delays;
probeAxis   = dataStruct.cmprobe;

Ndelays = length(delays);
Npixels = length(probeAxis);

switch dataStruct.rawcorr
    case 'CORRECTED'
        Z = dataStruct.corrdata;
    case 'RAW'
        Z = dataStruct.rawsignal;
end

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


mode        = questdlg('Select Chirp Correction Method','Chirp Correction Method','Automatic','Manual','Cancel','Automatic');

switch mode
    case 'Cancel'
        return
    case 'Manual'
        % Make a plot of the data
        maxPl = max(abs(Z(:)));
        minPl = -maxPl;
        ctrs  = linspace(minPl,maxPl,40);
        
        fh = figure(1);
        clf(fh);
        fh.Color = 'w';
        ax = axes('parent',fh);
        
        contourf(ax,probeAxis,delays,Z,ctrs,'EdgeColor','flat');
        
        colormap(darkb2r(minPl,maxPl,40,2));
        caxis(ax,[minPl,maxPl]);
        colorbar(ax);
        
        ylim([-5,5]); % Set this preliminary zoom level in case the data contains longer delays
        
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
        a=0;
        delete(fh);
    case 'Automatic'
        %%% Fit a Gaussian response and take the WL-dependent t0 parameter to do the chirp fit
        % 1D Gaussian with 1st and 2nd derivative
        % Parameters: 1=t0  2=FWHM  3=Amplitude  4=Offset  5=Amp 1st deriv  6=Amp 2nd deriv
        G0 = @(p,t) exp(-log(2).*((t-p(1))./p(2)).^2);
        G1 = @(p,t) -2.*log(2)./(p(2).^2).*(t-p(1)).*G0(p,t);
        G2 = @(p,t) 2.*log(2).*(-p(2).^2 + 2.*log(2).*(t-p(1)).^2)./(p(2).^4).*G0(p,t);

        ArtFit = @(p,t) p(3) + p(4).*G0(p,t) + p(5).*G1(p,t) + p(6).*G2(p,t);

        P0  = [0    0.1   10    0.1    0.1  0.1  ];
        LB  = [-5   0     -50   -50  -50  -50  ];
        UB  = [5    1     50    50   50   50   ];

        % Prepare an array where we will store the results
        Pfit = zeros(NfitPix,length(P0));
        Dfit = zeros(Ndelays,NfitPix);

        % Define fit options
        options     = optimoptions(@lsqcurvefit,...
                        'FunctionTolerance',5e-10,...
                        'MaxIterations',1e4,...
                        'MaxFunctionEvaluations',1e4,...
                        'steptolerance',5e-10,...
                        'display','off');

        % Do the fits
        wb = waitbar(0);            
        for j=1:NfitPix
            i=fitPixels(j);
            [P0(3),idM] = max(abs(Z(:,i)));
            P0(1)   = delays(idM);
            UB(1)   = P0(1) + 0.2;
            LB(1)   = P0(1) - 0.2;

            UB(3:6)   = 10*P0(3);
            LB(3:6)   = -10*P0(3);
            Pfit(j,:) = lsqcurvefit(ArtFit,P0,delays,Z(:,i),LB,UB,options);
            Dfit(:,j) = ArtFit(Pfit(j,:),delays);
            waitbar(j/NfitPix,wb,['Fitting chirp correction... (' num2str(j) ' of ' num2str(NfitPix) ')'])
        end
        delete(wb);
end

%% Fit Chirp vs Lambda

% Fit dispersion using Cauchy's Equation (truncated on the 3rd/4th term)
% Parameters: 1=a 2=b 3=c; (4=d)  L = Wavelength in nm
chirpFun = @(p,L) p(1) + 1e8.*p(2)./(L.^2) + 1e12.*p(3)./(L.^4);% + 1e17.*p(4)./(L.^6);

C0 = [3     1       1   ];%    1    ];
UB = [5     1e3     1e3  ];%    1e7  ];
LB = [-5    -1e3    -1e3 ];%    -1e7 ];

options     = optimoptions(@lsqcurvefit,...
                'FunctionTolerance',1e-10,...
                'MaxIterations',5e4,...
                'MaxFunctionEvaluations',5e4,...
                'steptolerance',1e-8,...
                'display','off');
            
Cfit = lsqcurvefit(chirpFun,C0,fitWL,Pfit(:,1),LB,UB,options);

% Plot Everything
%%% Plot fitted chirp data
WHAT = Z;
maxPl = max(abs(WHAT(:)));
minPl = -maxPl;
ctrs  = linspace(minPl,maxPl,40);

fh = figure(1);
clf(fh);
fh.Color = 'w';

ax = axes('parent',fh);

contourf(ax,probeAxis,delays,WHAT,ctrs,'EdgeColor','flat');
colormap(darkb2r(minPl,maxPl,40,2));
caxis([minPl,maxPl]);
colorbar;

hold(ax,'on');
plot(ax,fitWL,Pfit(:,1),'o','Color',0.8*[0 1 1])
plot(ax,probeAxis,chirpFun(Cfit,probeAxis),'y','LineWidth',3);
hold(ax,'off');

ylim(ax,[min(Pfit(:,1))-0.25 max(Pfit(:,1))+0.25])
title(ax,'Fitted Dispersion Curve','FontWeight','bold');
ylabel(ax,'Delay (ps)','FontWeight','bold');
xlabel(ax,'Wavelength (nm)','FontWeight','bold');
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

res = 1000*(Pfit(:,1)-chirpFun(Cfit,fitWL));

plot(ax2,fitWL,res,'xk')
yline(ax2,0);
ax.XTickLabel=[];
xlabel(ax2,'Wavelength (nm)','FontWeight','bold');
ylabel(ax2,'{\Delta}{t_{0}}^{res} (fs)','FontWeight','bold');
ylim(ax2,[-50,50]);

box(ax2,'on');
ax.FontSize = 16;
ax2.FontSize = 16;
ax.TickLength = [0.02 0.02];
ax2.TickLength = [0.02 0.02];
linkaxes([ax,ax2],'x');

switch mode
    case 'Automatic'
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

%% Save the chirp correction .mat file
path_save = uigetdir(rootfolder,'Save Chirp Correction Parameters to...');
fn_save = ['chirpCorr_' datestr(datetime('now'),'yyyymmdd-HHMMSS') '.mat'];
try 
    save([path_save filesep fn_save],'chirpFun','Cfit','fitPixels','fitWL','Pfit');
    helpdlg(["Chirp correction file successfully saved to: "; [path_save filesep fn_save]],'Chirp Correction Saved!');
catch ME
    errordlg('Error saving chirp correction file','Error saving file');
end

%% Copy variables to output structure
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