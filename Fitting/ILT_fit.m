function [Tau,LTspec] = ILT_fit(time,probeWL,Zdata,Noise)
% Inverse Laplace Transform (Regularised by maximum entropy method)
% Preliminary version (0.9b)
% 31.07.2020 / Ricardo Fernández-Terán
%

plotWL      = 20;
test        = 0;

%% Hardcoded Settings
n10         = 15;                       % No. time constants per decade
tau_min     = 0.1;                      % Shortest time constant
tau_max     = 1e5;                      % Longest time constant
% weight      = mean(abs(Noise(:)));      % Weight of the entropy contribution
weight = 1e-3;
% weight = 1.4709;


%% Generate dummy data
if test==1
    Ntaus       = round(log(tau_max/tau_min)/log(10)*n10)+1;
    Tau         = logspace(log10(tau_min),log10(tau_max),Ntaus);
    time = logspace(log10(tau_min),log10(tau_max),100);
    
    wl = linspace(350,750,32);

    St = 0*exp(-((Tau-1).^2)/50) - 0.5*exp(-((Tau-80).^2)/250) + 0.8.*exp(-((Tau-800).^2)/2500) - 0.6.*exp(-((Tau-2500).^2)/7500);
    Wt = -1.5.*exp(-(wl-420).^2/2500) + 0.75.*exp(-(wl-620).^2/10000) + 0.5.*exp(-(wl-580).^2/10000);
    Wt = Wt';
    Dt = Wt*St;

    kin = St*exp(-time'./Tau)';
    DAbs = Wt*kin;

    probeWL = wl';
    Zdata   = DAbs';
    Noise   = 0.0001/100.*max(abs(Zdata(:))).*rand(size(Zdata,1),size(Zdata,2));
    Zdata   = Zdata + Noise;
end


%% Initialise Variables
Ntaus       = round(log(tau_max/tau_min)/log(10)*n10)+1;
Tau         = logspace(log10(tau_min),log10(tau_max),Ntaus);
Ntime       = length(time);
Nprobe      = length(probeWL);
RMS_n       = rms(Noise(:));

% weight      = RMS_n*weight;

if size(time,2) == 1
    time        = time';
end
% plot(time,kin);
% ax=gca;
% ax.XScale = 'log';
% 
% contourf(time,wl,DAbs,50);
% ax.XScale = 'log';

%% Define Fit Options
options     = optimoptions(@fminunc,'PlotFcn',@optimplotfval,'MaxIterations',50,'MaxFunctionEvaluations',Ntaus*Nprobe*200,'StepTolerance',1e-10,'FiniteDifferenceType','forward','UseParallel',true);

%% Write fit function
m           = (max(Zdata(:)) - min(Zdata(:)))/Ntaus/100;
StartParam  = zeros(Ntaus+1,Nprobe);

[LTspec,fval,exitflag,output,grad] = fminunc(@(x)FitFunc(x,Zdata,time,Tau,m,weight),StartParam,options);

%%
options     = optimoptions(@fminunc,'PlotFcn',@optimplotfval,'MaxIterations',100,'MaxFunctionEvaluations',Ntaus*Nprobe*200*10,'StepTolerance',1e-10,'FiniteDifferenceType','central','TypicalX',LTspec,'UseParallel',true);
[LTspec,fval,exitflag,output,grad] = fminunc(@(x)FitFunc(x,Zdata,time,Tau,m,weight),LTspec,options);


%% Reconstruct data
Zfit    = exp(-time'./Tau)*LTspec(2:end,:) + LTspec(1,:);
Res     = Zfit-Zdata;
RMS_f   = sum(Res(:).^2)/Ntaus/Nprobe;

disp(['RMS data noise: ' num2str(RMS_n)])
disp(['RMS residuals:  ' num2str(RMS_f)])

%% Other plots
DT      = sqrt(sum(LTspec.^2,2));     DT = DT/max(DT(:));
% DT_th   = sqrt(sum(Dt.^2,1));         DT_th = DT_th/max(DT_th(:));

fh3 = figure(3); clf(fh3);
ax3 = axes('parent',fh3);
    plot(Tau,DT(2:end));
%     hold(ax3,'on');
%     plot(Tau,DT_th);
ax3.XScale = 'log';

Ncontours = 50;
%
fh1 = figure(1); clf(fh1);
ax1 = axes('parent',fh1);
    contourf(ax1,Tau,probeWL,LTspec(2:end,:)',Ncontours); colorbar
%     contourf(ax1,time,probeWL,Zfit',Ncontours); colorbar
ax1.XScale = 'log';    

fh2 = figure(2); clf(fh2);
ax2 = axes('parent',fh2);
    contourf(ax2,Tau,probeWL,Dt,Ncontours); colorbar
%     contourf(ax2,time,probeWL,(Zfit-Zdata)',Ncontours,'EdgeColor','flat'); colorbar
%     contourf(ax2,time,probeWL,Zdata',Ncontours); colorbar
    title(ax2,'Real / residual');
ax2.XScale = 'log';

%% Combined plot
% Kinetics
plotWL_id = unique(findClosestId2Val(probeWL,plotWL));
Nplots    = length(plotWL_id);

for i=1:Nplots
    % Create figure
    fh(i) = figure(3+i); clf(fh(i));
    ax(i) = axes('parent',fh(i));
    hold(ax(i),'on');
%     yyaxis(ax(i),'left');
    % Plot raw data
    plot(ax(i),time,Zdata(:,plotWL_id(i)),'ok');
    % Plot fitted kinetics
    plot(ax(i),time,Zfit(:,plotWL_id(i)),'-r','LineWidth',2);
    % Plot fitted lifetime spectrum
%     yyaxis(ax(i),'right');
    plot(ax(i),Tau,LTspec(2:end,plotWL_id(i)),'-b','LineWidth',2);
%     yyaxis(ax(i),'left');
    % Axis scale
    ax(i).XScale = 'log';
    hold(ax(i),'off');
    yline(ax(i),0,'HandleVisibility','off');
end

Tau= Tau';

%%
function out = FitFunc(x,Z,t,T,m,Lambda)
Zfit = exp(-t'./T)*x(2:end,:) + x(1,:);

RMSD = norm(Z-Zfit)./length(t);
S    = sum(sqrt(x.^2+4*m.^2) - x.*log((sqrt(x.^2+4*m.^2)+x)./(2.*m)) - 2.*m,'all');

out  = RMSD./Lambda - S;
end

end