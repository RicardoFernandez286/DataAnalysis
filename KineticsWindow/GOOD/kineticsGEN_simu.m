function Ct = kineticsGEN_simu(t,k,C0,ODEfunc,MakePlot,NameArray)
% Description:  This function returns a concentration matrix (Ct) of size nt x N for a system with a given time vector (t)
%               initial concentrations (C0), obeying an ODE system (ODEfunc) with rate constants in an array (k).
% 
% Usage:        Ct = kineticsGEN_simu(t,K,C0,ODEfunc,MakePlot,NameArray)
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

%% Numerically solve the differential equation
[~,Ct] = ode23tb(@(t,C0) ODEfunc(C0,k), t, C0, opt);  % Solve the ODE system and store the solution in the "Sol" matrix

%% Plot results
if MakePlot==true
    % Get input sizes
    N   = length(C0);

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
    
    % Show X axis in logarithmic scale and set proper limits
    ax.XScale   = 'log';
    xlim(ax,[min(t),max(t)]);

    % Resize Figure and Axes
    ax.Units    = 'pixels';
    fh.Units    = 'pixels';
    fh.Position(3:4) = [1072 420];
    ax.Position = [70 60 960 340];
    ax.Units    = 'normalized';
    fh.Units    = 'normalized';
end