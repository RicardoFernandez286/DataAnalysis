function Sv = EnT_kinetics(P, t)
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
x0 = [P(5) P(6) P(7) P(8)]; 

% Numerically solve the differential equation
[~,Sv] = ode23s(@DiffEq, t, x0);

% Define the system of differential equations to solve
function dS = DiffEq(t, x)
  n_eqs = 4;
%   dS = zeros(n_eqs,1);
    A     = - x(1)./P(1) - x(1)./P(3) + x(2)./P(4);
    AB    = - x(2)./P(2) - x(2)./P(4) + x(1)./P(3);
    B     = - x(3)./P(1) - x(3)./P(4) + x(4)./P(3);
    BA    = - x(4)./P(2) - x(4)./P(3) + x(3)./P(4);
  dS = [A;AB;B;BA];
%   dS = [A;AB];
end

end