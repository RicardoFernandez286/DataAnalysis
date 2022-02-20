function Ct = kineticsKmat_simu(t,K,C0,MakePlot,NameArray)
% Description:  This function returns a concentration matrix (Ct) of size nt x N for a system with a given time vector (t)
%               initial concentrations (C0) and rate matrix (K).
% 
% Usage:        Ct = kineticsKmat_simu(t,K,C0,MakePlot,NameArray)
%
% Inputs:       t = time axis (nt x 1);  K = rate matrix (must be square, NxN);  C0 = concentration vector (N x 1)
%               MakePlot = true or false (boolean) to decide whether a plot of the results is wanted or not.
%               NameArray = a string array with the names of the species to be plotted (optional).
%
% Outputs:      Ct = concentration matrix (nt x N)
%
% Tested and implemented in MATLAB R2021b
% v1.0

%% Build Concentration Matrix
EigenVal = 1;   % Set to zero in order to use the matrix exponential calculation,
                % otherwise use the eigenvalue/eigenvector method of Berberan-Santos and Martinho.

t   = t(:)'; % Ensure t is always a ROW vector
C0  = C0(:); % Ensure C0 is always a COLUMN vector

% Get input sizes
nt  = length(t);
N   = length(C0);

if EigenVal == 1
    [V,Lambda] = eig(K);
    Lambda=diag(Lambda);
    a = V\C0;
    ExpVec = exp(Lambda*t);
    Ct = (V*(a.*ExpVec))';
else
    % Preallocate for speed
    Ct  = zeros(nt,N);
    % Calculate the matrix exponentials for every time point
    for i=1:length(t)
        Ct(i,:) = expm(K*t(i))*C0;
    end
end
%% Plot results
if MakePlot==true
    % Create a new figure and axes
    fh      = figure(1);
    clf(fh)
    fh.Color= 'w';
    ax      = axes('Parent',fh);
    cmap    = brighten(turbo(N+1),-0.25);
    colororder(ax,cmap);
    SpeciesNames = cell(N,1);
    % Plot the curves
    hold(ax,'on')
    for j=1:N
        if isempty(NameArray)
            SpeciesNames{j} = ['Species ' num2str(j)];
        else
            SpeciesNames = NameArray;
        end

        plot(t,Ct(:,j),'linewidth',2,'DisplayName',SpeciesNames{j})
    end
    hold(ax,'off')
    
    % Add legend and box
    legend(ax,'show','Location','best','interpreter','latex','box','off');
    box(ax,'on');
    
    % Add axes labels, set font and tick size
    xlabel(ax,'Time','FontWeight','bold');
    ylabel(ax,'Concentration','FontWeight','bold');
    ax.FontSize     = 14;
    ax.TickLength   = [0.02 0.02];
 
    % Resize Figure and Axes
    ax.Units    = 'pixels';
    fh.Units    = 'pixels';
    fh.Position(3:4) = [1072 420];
    ax.Position = [70 60 960 340];
    ax.Units    = 'normalized';
    fh.Units    = 'normalized';
end