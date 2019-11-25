maxdelay=30;

HeNe    = 2.11079*10^-3;          % HeNe period (fs)
c       = 2.99792458e-2;    % Speed of light in cm/ps
t       = linspace(-maxdelay,maxdelay,2*maxdelay/HeNe); % ps
w1      = linspace(0,1/(c*HeNe),length(t));

omega = 2100; %cm^-1
FWHM  = 0.150; %ps

sigma = FWHM/sqrt(2*log(2));

d=-0.5; %ps

pulse1=exp(-0.5*(t/sigma).^2).*cos(2*pi*omega*c*(t));
pulse2=0.1*exp(-0.5*((t-d)/sigma).^2).*cos(2*pi*omega*c*(t-d));

figure(1)
plot(t,pulse1,'b')
hold on
plot(t,pulse2,'r')
hold off

% FFT_p1 =abs(fft(pulse1));
% FFT_p2 =abs(fft(pulse2));
FFT_sum=(abs(fft(pulse1+pulse2))).^2;

% figure(2)
% % hold on
% plot(w1,FFT_sum);
% xlim([1900,2250])
% % hold on
% plot(w1,FFT_p1);
% hold off

%% Pepita
figure(3)
hold off
pp_signal = log(abs(fft(pulse1+0.1*pulse1))./abs(fft(pulse1+pulse2)));
plot(abs(fft(pp_signal(3200:4300)-mean(pp_signal(3200:4300)))))
% plot(pp_signal(3200:4300)-mean(pp_signal(3200:4300)))
% xlim([3200,4400])
xlim([0 100]);