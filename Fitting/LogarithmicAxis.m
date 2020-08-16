tmin    = 0.2;
tmax    = 20000;
dT      = 0.2;

Nt      = 43;

T       = logspace(log10(round(tmin/dT)),log10(round(tmax/dT)),Nt)';
T       = unique(round(T))*dT;
T       = [0; T];