function FitProbeSpectrum(CAL_data)

%% Get the data and normalise it first
I1Norm = CAL_data.Air{1}./max(CAL_data.Air{1});
I2Norm = CAL_data.Air{2}./max(CAL_data.Air{2});
cm = CAL_data.cm;

%% Fit Gaussians
[~,pos1] = max(I1Norm);
[~,pos2] = max(I2Norm);

p01 = [1 cm(pos1) 250];
p02 = [1 cm(pos2) 250];

LB  = [0 0 1];
UB  = [5 4000 1000];

opt            = optimset(@lsqcurvefit);
opt.TolFun     = 1e-3;
opt.MaxIter    = 150;
opt.MaxFunEvals= 500;
opt.TolX       = 1e-3;
opt.FinDiffType='central';
opt.Display    = 'off';

G1D = @(P,X) P(1).*exp(-4*log(2)*((X-P(2))./P(3)).^2);

pFit1 = robustlsqcurvefit(G1D,p01,cm(:,1),I1Norm,LB,UB,'bisquare',opt);
pFit2 = robustlsqcurvefit(G1D,p02,cm(:,2),I2Norm,LB,UB,'bisquare',opt);

cm_fit = linspace(min(cm(:))-500,max(cm(:))+500,1000);
g1_plt = G1D(pFit1,cm_fit);
g2_plt = G1D(pFit2,cm_fit);

I1Norm = I1Norm./pFit1(1);
I2Norm = I2Norm./pFit2(1);

g1_plt = g1_plt./pFit1(1);
g2_plt = g2_plt./pFit2(1);

meanW0 = mean([pFit1(2),pFit2(2)]);
minX = meanW0 - 1.5*mean([pFit1(3),pFit2(3)]);
maxX = meanW0 + 1.5*mean([pFit1(3),pFit2(3)]);

%% Plot Results
fhS = figure(3);
fhS.Color = 'w';
fhS.Position(3:4) = [600 350];
clf(fhS);
axS = axes('parent',fhS);

hold(axS,'on');
    plot(axS,cm(:,1),I1Norm,'o','Color',[0 0.5 0],'MarkerSize',4,'DisplayName','Det 1');    
    plot(axS,cm(:,2),I2Norm,'o','Color',[0  0  1],'MarkerSize',4,'DisplayName','Det 2');

    plot(axS,cm_fit,g1_plt,'Color',[0 0.5 0],'LineWidth',1,'DisplayName','Det 1 (Fit)')
    plot(axS,cm_fit,g2_plt,'Color',[0  0  1],'LineWidth',1,'DisplayName','Det 2 (Fit)')
hold(axS,'off');


xlabel(axS,'Wavenumbers (cm^{-1})','FontWeight','bold');
ylabel(axS,'Norm. Intensity','FontWeight','bold');

axS.Box = 'on';
axS.FontSize = 14;

yline(axS,0.5,'--','HandleVisibility','off');
xline(axS,meanW0,'--','HandleVisibility','off');
yline(axS,0,'-','HandleVisibility','off');

xlim(axS,[minX maxX])
S = sprintf('FWHM = %.4g / %.4g cm^{-1};  w0 = %.4g / %.4g cm^{-1}',pFit1(3),pFit2(3),pFit1(2),pFit2(2));

% fprintf('\nGaussian Fit Info (D%i):\n FWHM = %.4g cm-1\n w0   = %.4g cm-1\n\n',1,pFit1(3),pFit1(2))
% fprintf('Gaussian Fit Info (D%i):\n FWHM = %.4g cm-1\n w0   = %.5g cm-1\n',2,pFit2(3),pFit2(2))

title(axS,S,'FontSize',14)

end
