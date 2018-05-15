
%% Sellmeier eq. calc. of the refractive index at the given wavelength

syms L % L = Wavelength (in ?m)

% Sellmeier equation for CaF2
n_CaF2(L)=sqrt(1+0.5675888./(1-(0.050263605./L).^2)+0.4710914./(1-(0.1003909./L).^2)+3.8484723./(1-(34.649040./L).^2));

% Refractive index of Si @ 5 ?m
n_Si=3.4195;

% Sellmeier equation for Air
% n_Air(L)=;

% Sellmeier equation for MeOH
n_MeOH(L)=1.294611+12706.403E-6.*L.^-2;

% Water at 5 ?m
n_H2O = 1.3145;

%% Select the refractive indexes:
% Wavelength (nm)
%     WL = 632.8/1000;

% Wavenumber (cm-1)
WN = 2000;

WL = 10000./WN; % WL in ?m
% Prism material
    n_1=n_Si;
% Solvent layer on top of the prism
    n_2=n_MeOH(WL);

%% Based on the two selected refractive indexes, calculate the TR angle:
    Theta_C=double(asind(n_2/n_1))