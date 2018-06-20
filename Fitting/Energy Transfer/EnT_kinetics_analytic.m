function [A_func,aB_func,bA_func,B_func] = EnT_kinetics_analytic(time,fw,bw)
% EnT_kinetics codes the system of differential equations describing energy transfer:

kB      = 8.6173303*10^-5; % eV K^-1
Temp    = 298; % K
dW      = 46; % cm-1
ev2cm   = 8065.54;
factor  = exp(-dW./((kB*Temp*ev2cm)));

syms A(t) aB(t) bA(t) B(t)
syms tau1 tau2 tauFW tauBW
syms C

k1 = 1./tau1;
k2 = 1./tau2;
kFW = 1./tauFW;
% kBW = kFW./factor;
kBW = 1./tauBW;

% FIRST DECAY
eqns_fw(1)     = diff(A,t)     == -k1*A  - kFW*A  + kBW*aB;
eqns_fw(2)     = diff(aB,t)    == -k1*aB - kBW*aB + kFW*A;
eqns_bw(1)     = diff(bA,t)    == -k1*bA - kFW*bA + kBW*B;
eqns_bw(2)     = diff(B,t)     == -k1*B  - kBW*B  + kFW*bA;

cond_fw(1)     = A(time(1))  == fw(1,1);
cond_fw(2)     = aB(time(1)) == 0;
cond_bw(1)     = bA(time(1)) == 0;
cond_bw(2)     = B(time(1))  == bw(1,1);

solFW = dsolve(eqns_fw,cond_fw);
solBW = dsolve(eqns_bw,cond_bw);

% SECOND DECAY
eqns_fw2(1)     = diff(A,t)     == -k2*A  - kFW*A  + kBW*aB;
eqns_fw2(2)     = diff(aB,t)    == -k2*aB - kBW*aB + kFW*A;
eqns_bw2(1)     = diff(bA,t)    == -k2*bA - kFW*bA + kBW*B;
eqns_bw2(2)     = diff(B,t)     == -k2*B  - kBW*B  + kFW*bA;

cond_fw2(1)     = A(time(1))  == fw(1,1);
cond_fw2(2)     = aB(time(1)) == 0;
cond_bw2(1)     = bA(time(1)) == 0;
cond_bw2(2)     = B(time(1))  == bw(1,1);

solFW2 = dsolve(eqns_fw2,cond_fw2);
solBW2 = dsolve(eqns_bw2,cond_bw2);

% FUNCTIONS
A_func(tau1,tau2,tauFW,tauBW,C,t)   = C.*solFW.A   + (1-C).*solFW2.A;
aB_func(tau1,tau2,tauFW,tauBW,C,t)  = C.*solFW.aB  + (1-C).*solFW2.aB;
B_func(tau1,tau2,tauFW,tauBW,C,t)   = C.*solBW.B   + (1-C).*solBW2.B;
bA_func(tau1,tau2,tauFW,tauBW,C,t)  = C.*solBW.bA  + (1-C).*solBW2.bA;

