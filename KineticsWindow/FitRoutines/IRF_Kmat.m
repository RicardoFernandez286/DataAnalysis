function [CConc,tIRF,simIRF,NConc] = IRF_Kmat(t,K,C0,w,t0)
% Description:  This function returns an IRF-convolved concentration matrix (CConc) of size nt x N for a system with a 
%               given time vector (t) and initial concentrations (C0), obeying a 1st order kinetic system with 
%               a given K matrix (K). The Gaussia
% 
% Usage:        [CConc,tIRF,simIRF,NConc] = IRF_Kmat(t,K,C0,w,t0)
%
% Inputs:       t = time axis (nt x 1);  k = array of rate constants;  C0 = initial concentration vector (N x 1).
%               OFEfunc = a function handle to a system of differential equations, taking concentration and time inputs.
%               MakePlot = true or false (boolean) to decide whether a plot of the results is wanted or not.
%               NameArray = a string array with the names of the species to be plotted (optional).
%
% Outputs:      Ct = concentration matrix (nt x N)
%
% Tested and implemented in MATLAB R2021b
% v1.0
% 

% Gaussian IRF function (area-normalised Gaussian function)
IRFfnc      = @(T,T0,W) (2./W)*sqrt(log(2)/pi)*exp(-4*log(2).*((T-T0)./W).^2)';

t           = t(:)'; % Ensure t is always a ROW vector
C0          = C0(:); % Ensure C0 is always a COLUMN vector

% Eigenvalue method of Berberan-Santos and Martinho for a 1st order kinetic system
[V,Lambda]  = eig(K);
Lambda      = diag(Lambda);
a           = V\C0;

% Exponential function convolved with a Gaussian IRF of given FWHM (w) and t0.
s           = w/(2*sqrt(2*log(2))); % Convert FWHM to sigma, since it gives simpler equations
ExpVec_GC   = 0.5.*exp(Lambda.*(t - t0 + 0.5*Lambda*s.^2)).*(1+erf((t-t0+Lambda.*s.^2)./(s.*sqrt(2))));

% Convolved response, as the product of ExpVec and the corresponding pre-exp. factors from C0 and the eigenvectors of K 
CConc       = (V*(a.*ExpVec_GC))';

% Non-convolved response 
ExpVec_Nor  = exp(Lambda.*(t-t0));
NConc       = (V*(a.*ExpVec_Nor))';
NConc       = NConc((t-t0)>=0,:);

% Generate the gaussian for the convolution
endConv     = 100*w;
step        = min(diff(t));
timeConv    = (min(t)-rem(min(t),step):step:endConv);
tIRF        = (-ceil(length(timeConv)*0.4)*step:step:ceil(length(timeConv)*0.4)*step)';
simIRF      = IRFfnc(tIRF,t0,w);
