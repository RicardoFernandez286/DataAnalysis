function grTot = calcGR(trajectData_3D,Nbins,BoxSizes,Nmolec,PBC)

% Initialise certain parameters
Nframes = size(trajectData_3D,3);   % Number of MD frames
rc      = min(BoxSizes)/2;          % Cutoff radius  
dr      = rc/(Nbins+1);             % Bin size
g_x     = (1:Nbins)'*dr;
rho     = Nmolec/prod(BoxSizes);    % Particle density

% Get only the XY coordinate vectors
pos     = trajectData_3D(:,1:2,:);

%% Accumulate counts
g       = zeros(Nbins+1,Nframes);
for f=1:Nframes
% Count
    for n1=1:Nmolec
        rij = (pos(n1,:,f) - pos(:,:,f));
        rij = rij - round(rij./BoxSizes).*BoxSizes*PBC;
        dij = sqrt(sum(rij.^2,2));
        dij = dij(dij<rc);
        dij = dij(dij>0);
        idx = ceil(dij/dr);
        for ij=1:length(dij)
            g(idx(ij),f) = g(idx(ij),f) + 1;
        end
    end
    % Normalise
    for n=1:Nbins
        g(n,f)  = g(n,f)/Nmolec;        % Divide by Nmolec
        dA      = 2*pi*(dr*n)*dr;       % Area of a disc (2D)
%       dV      = 4*pi*(dr*n)^2*dr;     % Volume of a spherical shell (3D)
        g(n,f)  = g(n,f)/dA;            % g = Local density
        g(n,f)  = g(n,f)/rho;           % g = RDF
    end
end
g       = mean(g(1:end-1,:),2);

grTot   = [g_x, g];