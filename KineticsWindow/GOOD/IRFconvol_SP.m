function [CConc,tIRF,simIRF] = IRFconvol_SP(IRFfnc,t,t0,w,ConcFnc)
% Description:  This function returns a concentration matrix (CConc) of size nt x N for a system with a concentration function
%               (ConcFnc) that is evaluated in time (t) and convolved with a given IRF (IRFfnc) of time zero and FWHM given by
%               t0 and w, respectively.
% 
% ** Special version to take into account a positive concentration of the initial reagent at t<t0 **
% 
% Usage:        [CConc,tIRF,simIRF] = IRFconvol(IRFfnc,t,t0,w,ConcFnc)
%
% Inputs:       IRFfnc = a function handle to the IRF, IRF(t,t0,w); where t = time axis (nt x 1), t0 = time zero, w = FWHM
%               ConcFnc = a function handle to the concentration function. It should only depend on the time axis [i.e. C(t)]
%
% Outputs:      CConc = IRF-convolved time-dependent concentration matrix
%               tIRF and simIRF = time and IRF evaluated at the given times, in case it is to be plotted with the data.
%
% Tested and implemented in MATLAB R2021b
% v1.0

%%% Initial definitions
t           = t(:)';    % Ensure t is always a ROW vector
endConv     = 5000*w;   % Changed so it works better with the kinetic data for the chameleon reaction
step        = min(diff(t));

timeConv    = (min(t)-rem(min(t),step):step:endConv);   % time range where the convolution is done
timeSim     = [timeConv t(t > endConv)];                % time axis for the non-convolved response 

% Heaviside step function, usual definition with H(0) = 1/2
Hfunc       = @(t) (1+sign(t))/2;

%%% Calculate non-convolved responses
% Concentrations are multiplied by the Heaviside step function since they must be 0 at tSim<0
C_H         = Hfunc(timeSim)'.*ConcFnc(timeSim);

%!! Positive concentration of the initial reagent at t<t0
C_H(timeSim<=0,1) = 1;

%%% Convolution with Gaussian IRF
% Generate the Gaussian IRF for the convolution
tIRF        = (-ceil(length(timeConv)*0.4)*step:step:ceil(length(timeConv)*0.4)*step)';
simIRF      = IRFfnc(tIRF,t0,w).*step;

% Do the convolution
conv_full   = conv2(simIRF,1,C_H(timeSim <= endConv,:),'same');

% Shift the time of the convolved spectra because of the IRF
timeConv    = timeConv-t0;

% Time after which the non-convolved response is used instead
% Needs to be smaller than endConv to avoid getting 
% an artificial decay due to truncation in the numerical convolution
cutoff      = endConv-10*w; 

% Add the the two parts together with the correct time axes
time_full   = [timeConv(timeConv<cutoff), timeSim(timeSim>=cutoff)];
conv_full   = [conv_full(timeConv<cutoff,:); C_H(timeSim>=cutoff,:)];

% Finally, interpolate the concentration on the original time axis accounting for the shift of the IRF
CConc       = interp1(time_full+t0,conv_full,t,'spline'); % 'spline' gave best results

% Ensure that CConc is a COLUMN vector if we have only 1 species
if size(CConc,1) == 1
    CConc = CConc';
end

%%!! Positive concentration of the initial reagent at t<t0
CConc(t<(t0-w),1) = 1;

end