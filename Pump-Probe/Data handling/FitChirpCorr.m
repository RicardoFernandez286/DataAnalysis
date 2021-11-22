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
probeRanges = [0 1000];

%% Select data
ProbeIdx    = findClosestId2Val(probeAxis,probeRanges);
fitPixels   = [];

for j=1:round(length(probeRanges)/2)
    fitPixels = [fitPixels ProbeIdx(2*(j-1)+1):1:ProbeIdx(2*j)];
end
fitPixels   = unique(fitPixels);
fitWL       = probeAxis(fitPixels);

NfitPix     = length(fitPixels);

%% Fit a Gaussian response and take the WL-dependent t0 parameter to do the chirp fit

% 1D Gaussian with offset as fit function
% Parameters: 1=t0 2=FWHM 3=Amplitude 4=Offset
G0 = @(p,t) p(3).*exp(-log(2).*((t-p(1))./p(2)).^2) + p(4);
G1 = @(p,t) -2.*log(2)./p(2).*t.*G0(p,t);
G2 = @(p,t) -2.*log(2)./p(2).*G0(p,t) - G1(p,t).*2.*log(2)./p(2).*t;

% syms t t0 A s A0 A1 A2
% Gau1D_sym(t,t0,A,s) = A.*exp(-log(2).*((t-t0)./s).^2);
% Gau1D_art(t,t0,s,A,A0,A1,A2) = Gau1D_sym(t,t0,A,s) + A0 + A1*diff(Gau1D_sym(t,t0,A1,s),1) + A2*diff(Gau1D_sym(t,t0,A2,s),2);
% ArtFit = @(p,t) (Gau1D_art(t,p(1),p(2),p(3),p(4),p(5),p(6)));

ArtFit = @(p,t) G0(p,t) + G1(p,t) + G2(p,t);

% Initial Parameters and bounds
% P0  = [0    0.15  10    0   ];
% LB  = [-5   0     -50   -50 ];
% UB  = [5    1     50    50  ];

P0  = [0    0.15  10    10   10   10  ];
LB  = [-5   0     -50   -50  -50  -50];
UB  = [5    1     50    50   50   50 ];


% Prepare an array where we will store the results
Pfit = zeros(NfitPix,length(P0));

% Define fit options
options     = optimoptions(@lsqcurvefit,...
                'FunctionTolerance',5e-3,...
                'MaxIterations',5e2,...
                'MaxFunctionEvaluations',5e2,...
                'steptolerance',5e-3,...
                'display','off');

% Do the fits
wb = waitbar(0);            
for j=1:NfitPix
    i=fitPixels(j);
    [~,idM] = max(abs(Z(:,i)));
    P0(1)   = delays(idM);
    Pfit(j,:) = lsqcurvefit(ArtFit,P0,delays,Z(:,i),LB,UB,options);
    waitbar(j/NfitPix,wb,['Fitting chirp correction... (' num2str(j) ' of ' num2str(NfitPix) ')'])
end
delete(wb);

%% Fit Chirp vs Lambda

% Fit dispersion using Cauchy's Equation (truncated on the 3rd/4th term)
% Parameters: 1=a 2=b 3=c; (4=d)  L = Wavelength in nm
chirpFun = @(p,L) p(1) + 1e7.*p(2)./(L.^2) + 1e11.*p(3)./(L.^4);% + 1e16.*p(4)./(L.^6);

C0 = [3     1       -1];%      1    ];
UB = [5     1e3     1e3];%      1e7  ];
LB = [-5    -1e3    -1e3];%     -1e7 ];

options     = optimoptions(@lsqcurvefit,...
                'FunctionTolerance',1e-10,...
                'MaxIterations',5e4,...
                'MaxFunctionEvaluations',5e4,...
                'steptolerance',1e-8,...
                'display','off');
            
Cfit = lsqcurvefit(chirpFun,C0,fitWL,Pfit(:,1),LB,UB,options);

%%% Plot dispersion fit
fh = figure(3);
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

%%% Plot wavelength-dependent FWHM
fh = figure(2);
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

%% Save the chirp correction .mat file
path_save = uigetdir(rootfolder,'Save Chirp Correction Parameters to...');
fn_save = ['chirpCorr_' datestr(datetime,'YYYYMMDD-hhmmss') '.mat'];
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
chirpSt.meanIRF     = meanIRF;
