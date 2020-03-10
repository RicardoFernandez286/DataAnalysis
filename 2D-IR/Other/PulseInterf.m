minT_i = -0.5;
maxT_i = 5;
dT_i   = 0.1;

maxdelay=50;
Nloops = 20;

HeNe    = 2.11079*10^-3;                                % HeNe period (fs)
c       = 2.99792458e-2;                                % Speed of light in cm/ps
t       = linspace(-maxdelay,maxdelay,2*maxdelay/HeNe); % ps
w1      = linspace(0,1/(c*HeNe),length(t));             % cm-1

omega = 2100; %cm^-1
FWHM  = 0.150; %ps

sigma = FWHM/sqrt(2*log(2));

for N=1:Nloops
    if mod(N,2)==1
        minT = minT_i;
        maxT = maxT_i;
        dT   = dT_i;
    else
        minT = maxT_i;
        maxT = minT_i;
        dT   = -dT_i;
    end
        
for d=minT:dT:maxT
% d=0.5; %ps

pulse1=exp(-0.5*(t/sigma).^2).*cos(2*pi*omega*c*(t));
pulse2=1*exp(-0.5*((t-d)/sigma).^2).*cos(2*pi*omega*c*(t-d));

if N==1 && d==minT
    fh1 = figure(1);
    clf(fh1);
    ax1 = axes('parent',fh1);
end
plot(ax1,t,pulse1,'Color',[0.5 0.5 1])
hold(ax1,'on')
plot(ax1,t,pulse2,'Color',[1 0.5 0.5])
% plot(t,pulse1+pulse2,'Color',0.5*[1 1 1])
hold(ax1,'off')
xlim(ax1,sort([3*minT_i 1.5*min(maxT_i,maxdelay)]));
drawnow

% FFT_p1 =abs(fft(pulse1));
% FFT_p2 =abs(fft(pulse2));
FFT_sum=(abs(fft(pulse1+pulse2))).^2;
FFT_sum0=4*abs(fft(pulse1)).^2;

if N==1 && d==minT
    fh2 = figure(2);
    clf(fh2);
    ax2 = axes('parent',fh2);
end
plot(ax2,w1,FFT_sum0,'Color','k','LineWidth',2);
hold(ax2,'on')
 plot(ax2,w1,FFT_sum,'b');
 xlim(ax2,[2000,2200])
hold(ax2,'off')
drawnow
% plot(w1,FFT_p1);
% hold off



%% Pepita

if N==1 && d==minT
    fh3 = figure(3);
    clf(fh3);
    ax3 = axes('parent',fh3);
end

hold(ax3,'off')
pp_signal = log((abs(fft(pulse1+pulse2)).^2+4*abs(fft(pulse1)).^2)./abs(fft(pulse1)).^2);
plot(ax3,w1,pp_signal)
% plot(pp_signal(3200:4300)-mean(pp_signal(3200:4300)))
% xlim([3200,4400])
% xlim([0 100]);
xlim(ax3,[2000,2200])
drawnow

% pause(0.2)
end
end