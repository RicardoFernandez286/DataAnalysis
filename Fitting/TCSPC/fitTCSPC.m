%% Load data
datadir = 'F:\Cage_T_dependant_PrCN';

IRF_fn  = '20230523_IRF_440nm_excitation_50ns.asc';
data_fn = '20230523_RhodamineB_440nm_excitation_50ns.asc';

% Read files
IRF_xy  = readmatrix([datadir filesep IRF_fn],'FileType','text');
data_xy = readmatrix([datadir filesep data_fn],'FileType','text');

doFit = 1;
%% Preprocessing
BkgCts = 8;
tBkg   = 2;
FS = 14; % FontSize

% Remove NaNs (from the end)
IRF_xy  = IRF_xy(1:end-1,:);
data_xy = data_xy(1:end-1,:);

% Reject data with <1 cts
IRF_xy  = IRF_xy(data_xy(:,1)>1,:);
data_xy = data_xy(data_xy(:,1)>1,:);

% Extract time and counts
Y = data_xy(:,1);
IRF = IRF_xy(:,1);
t = 25.*linspace(0,1,length(Y));


% Set IRF to 0 if below BkgCts cts., "clean up IRF"
IRF(IRF<BkgCts) = 0;

% Normalise IRF to data
IRF = IRF.*max(Y)./max(IRF);

% Calculate normalised errors
sigma   = sqrt(Y);

% Calculate Background (from data)
bkg     = mean(Y(t<tBkg));

%% Create plotting axes
fh=figure(1);
fh.Color='w';
tl = tiledlayout(fh,5,1,'TileSpacing','compact');

ax(1)=nexttile; % Autocorrelation
ax(2)=nexttile; % Residuals
ax(3)=nexttile(tl,[3,1]); % Data+IRF+Fit


%% Plot the data
hold(ax(3),'on');
    plot(ax(3),t,IRF,'color',[0.5 0.5 0.5],'DisplayName','IRF')
    plot(ax(3),t,Y,'r','DisplayName','Data')

axis(ax(3),'tight');
%%
p0 = [-1, 150, 4.5];
LB = [-2 0 0];
UB = [+2 Inf Inf];

corrIFE = false;
nIFE    = 0;
Nexp    = 1;

if doFit == 1
    pOpt = lsqcurvefit(@(p,t) fitFun(p,t,IRF,bkg,Nexp,corrIFE,nIFE),p0,t,Y,LB,UB);
else
    pOpt = p0;
end

yFit = fitFun(pOpt,t,IRF,bkg,Nexp,corrIFE,nIFE);
res  = (Y - yFit)./sigma;
AC   = xcorr(res,res,length(t)-1);
% AC = AC((length(t)):end);
AC = flip(AC(1:length(t)));
AC = AC./max(AC);

%% Plot the rest
plot(ax(3),t,yFit,'k','DisplayName','Fit','LineWidth',3);
hold(ax(3),'off');

plot(ax(2),t,res,'b','LineWidth',1.5);
plot(ax(1),t,AC,'color',[0 0.75 0],'LineWidth',1.5);

% Beautify
for i=1:3
    box(ax(i),'on');
    ax(i).FontSize = FS;
end

ax(3).YScale = 'log';
linkaxes(ax,'x');

for i=1:2
    ax(i).XTickLabel = [];
end

ylabel(ax(1),'AC','FontWeight','bold'); yline(ax(1),0,'HandleVisibility','off');
ylabel(ax(2),'Res.','FontWeight','bold'); yline(ax(2),0,'HandleVisibility','off');
ylabel(ax(3),'Counts','FontWeight','bold');
xlabel(ax(3),'Delay (ns)','FontWeight','bold');

axis(ax(3),'tight');
ylim(ax(2),[-10 +10]);
ylim(ax(1),[-1 +1]);

% sprintf('Tau1 = %.6f ± %.6f ns',)
%%% EOF 
function yFit = fitFun(p,t,IRF,bkg,Nexp,corrIFE,nIFE)
    shift = p(1);
    a1 = p(2);
    tau1 = p(3);
    
    Npts = length(t);
    dT = mean(diff(t));
    IRF_int = interp1(t,IRF,t-shift,'linear',0);
    fitDec  = a1.*exp(-t./tau1);
    decC = conv(dT.*IRF_int,dT.*fitDec,'full');
%     yFit = conv2(IRF_int,1,dT.*fitDec);
    yFit = decC(1:Npts) + bkg;
    yFit = yFit(:);
end
