% function ILT_fit

%% Hardcoded settings

%% Initialise variables
t = logspace(log10(0.01),log10(42000),150);
wl = [350:1:750];

% Generate dummy data
St = -1.*exp(-((t-2).^2)/20) + exp(-((t-50).^2)/250) + 2.5.*exp(-((t-500).^2)/2500) - 1.5.*exp(-((t-1500).^2)/25000);
Wt = -1.5.*exp(-(wl-420).^2/2500) + 0.75.*exp(-(wl-620).^2/10000) + 0.5.*exp(-(wl-580).^2/10000);
Wt = Wt';
Dt = Wt*(-1.*exp(-((t-2).^2)/20)) + Wt*(exp(-((t-50).^2)/250)) + Wt*(2.5.*exp(-((t-500).^2)/2500)) + Wt*(- 1.5.*exp(-((t-1500).^2)/25000));

plot(wl,Wt);
ax=gca;
ax.XScale = 'log';

contourf(Dt,50);