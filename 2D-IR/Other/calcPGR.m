function [gx,gTot,gAA,gBB,gCC,gAB,gAC,gBC] = calcPGR(trajectData_3D,Nbins,BoxSizes,Nmolec,PBC,isoShift,varargin)

% Initialise certain parameters
Nframes = size(trajectData_3D,3);   % Number of MD frames
rc      = min(BoxSizes)/2;          % Cutoff radius  
dr      = rc/(Nbins+1);             % Bin size
gx      = (1:Nbins)'*dr;
% rho     = Nmolec/prod(BoxSizes);    % Particle density

% Get only the XY coordinate vectors
pos     = trajectData_3D(:,1:2,:);

% Classify particles according to their isotope shift ("1"=12, "2"=13; "3"=1318)
type    = zeros(size(trajectData_3D,1),size(trajectData_3D,3));
for i=1:3
    type(squeeze(trajectData_3D(:,6,:) == -isoShift(i))) = i;
end

% In the dilution with small molecules, the small molecules are classified as type 1, so they are shifted to type 3 instead
if ~isempty(varargin)
    r_small = varargin{1};
    type(squeeze(trajectData_3D(:,7,:) == r_small)) = 3;
    type(type==0) = NaN;
end
%% Accumulate counts
% msg = msgbox('Calculating pRDF. This may take a while...', 'Calc. in progress...','help','nonmodal');

msg     = waitbar(0,'Calculating pRDF. This may take a while...','WindowStyle','modal');
g       = zeros(Nbins+1,Nframes,9); % Nbins+1 x Nframes x 6 (AA,BB,CC,AB,AC,BC)

% The type is calculated from the product of the types of each molecule (i.e. 3x1 = 3 = BC)
%     AA = 1
%     BB = 4
%     CC = 9
%     AB = 2
%     AC = 3
%     BC = 6
molType = [1 4 9 2 3 6];
gpart   = zeros(Nbins,6);
for f=1:Nframes
    waitbar(f/Nframes,msg,'Calculating pRDF. This may take a while...','WindowStyle','modal');
    
    if length(unique(Nmolec)) > 1
        Nmolec_i = Nmolec(f);
    else
        Nmolec_i = unique(Nmolec);
    end
    rho_i    = Nmolec_i/prod(BoxSizes);
    for itype=1:6
        % Count
        for n1=1:Nmolec_i
            typIdx  = type(n1,f).*type(:,f);
            rij = (pos(n1,:,f) - pos(:,:,f));
            rij = rij - round(rij./BoxSizes).*BoxSizes*PBC;
            dij = sqrt(sum(rij.^2,2));
            dij = dij(typIdx==molType(itype));
            dij = dij(dij<rc);
            dij = dij(dij>0);
            idx = ceil(dij/dr);
            for ij=1:length(dij)
                g(idx(ij),f,itype) = g(idx(ij),f,itype) + 1;
            end
        end
        % Normalise
        for n=1:Nbins
            g(n,f,itype)  = g(n,f,itype)/Nmolec_i;        % Divide by Nmolec
            dA      = 2*pi*(dr*n)*dr;                   % Area of a disc (2D)
%             dV      = 4*pi*(dr*n)^2*dr;                 % Volume of a spherical shell (3D)
            g(n,f,itype)  = g(n,f,itype)/dA;            % g = Local density
            g(n,f,itype)  = g(n,f,itype)/rho_i;           % g = RDF
        end
        gpart(:,itype) = squeeze(sum(g(1:end-1,:,itype),2))/Nframes;
    end
end

gAA = gpart(:,1);
gBB = gpart(:,2);
gCC = gpart(:,3);
gAB = gpart(:,4);
gAC = gpart(:,5);
gBC = gpart(:,6);
gTot= sum(gpart,2);

close(msg);


