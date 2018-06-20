function S = EnT_kinetics(P, t, ~)
% EnT_kinetics codes the system of differential equations describing energy transfer:
%       dA/dt  = - k1*A  - k3*A  + k4*B
%       dAB/dt = - k2*AB - k4*AB + k3*A
%       dB/dt  = - k2*B  - k4*B  + k3*A
%       dBA/dt = - k1*BA - k3*BA + k4*B
%   with:
%       Variables:  x(1) = A,  x(2) = B
%       Parameters: P(1) = k1, P(2) = k2, P(3) = k3, P(4) = k4, P(5) = a0, P(6) = ab0, P(7) = b0,
%       P(8) = ba0.

% Define the initial conditions
x0 = [1 0 0 1];

% Numerically solve the differential equation
[~,S] = ode45(@DiffEq, t, x0);

% if nargin == 2
%     S = SX(:,1);
% elseif nargin == 3
%     S = SX;
% end

% Define the system of differential equations to solve
function dS = DiffEq(t, x)
%     n_eqs = 2;
%     dS = zeros(n_eqs,1);

    kB      = 8.6173303*10^-5; % eV K^-1
    Temp    = 298; % K
    dW      = -46; % cm-1
    ev2cm   = 8065.54;
    factor  = exp(-dW./((kB*Temp*ev2cm)));
    
    A       = (1-P(5)).*(- x(1)./P(1) - x(1)./P(3) + x(2)./P(4)) + P(5).*(- x(1)./P(2) - x(1)./P(3) + x(2)./P(4));
    aB      = (1-P(5)).*(- x(2)./P(1) - x(2)./P(4) + x(1)./P(3)) + P(5).*(- x(2)./P(2) - x(2)./P(4) + x(1)./P(3));
    bA      = (1-P(5)).*(- x(3)./P(1) - x(3)./P(3) + x(4)./P(4)) + P(5).*(- x(3)./P(2) - x(3)./P(3) + x(4)./P(4));
    B       = (1-P(5)).*(- x(4)./P(1) - x(4)./P(4) + x(3)./P(3)) + P(5).*(- x(4)./P(2) - x(4)./P(4) + x(3)./P(3));

    dS      = [A;aB;bA;B];
    
%     B     = - x(3)./P(1) - x(3)./P(4) + x(4)./P(3);
%     BA    = - x(4)./P(2) - x(4)./P(3) + x(3)./P(4);
%     dS = [A;aB;bA;B];
end

end