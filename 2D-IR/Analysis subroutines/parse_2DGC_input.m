function [PeaksFunction,Start_param,UB,LB,ParamPos] = parse_2DGC_input(fitparameters,Ndelays,varargin)
% Description: This function parses the fitparameters structure into a
% vector of starting parameters, upper bounds and lower bounds, together
% with the function to be used in the fit

%% Process the input structure

% Hardcoded settings
equal_SxSy  = 1;

if isempty(varargin)
    Omega   = {linspace(1900,2200,100);linspace(1900,2200,32)};
    ZData   = ones(length(Omega{1}),length(Omega{2}),Ndelays);
else
    Omega   = varargin{1};
    ZData   = varargin{2};
end

% Get the number of peaks (number of columns of the table) and sort them
% into Diagonal and Cross-peaks
Npeaks      = size(fitparameters,2);
isDiagonal  = contains(fitparameters(end,:),'Diag');

% Get the W1 and W3 axes
PumpAxis    = Omega{1};
ProbeAxis   = Omega{2};

% Get the starting parameters (as a cell)
X0_start    = cellfun(@str2num,fitparameters(1,:),'un',0); % Position along W1
Y0_start    = cellfun(@str2num,fitparameters(2,:),'un',0); % Anharmonicity/Position along W3 (see below)
SX_start    = cellfun(@str2num,fitparameters(3,:),'un',0); % Sigma_x
SY_start    = cellfun(@str2num,fitparameters(4,:),'un',0); % Sigma_y

% Convert the cells to double arrays
X0_start    = [X0_start{:}]';
Y0_start    = [Y0_start{:}]';
SX_start    = [SX_start{:}]';
SY_start    = [SY_start{:}]';
C_start     = fitparameters(5,:)';

% % The anharmonicities are smaller than 100 cm^-1. If a number >100 cm^-1 is
% % entered, then the user didn't mean anhmarmonicity but absolute position
% % of the peak along W3. This is stored in the vector "isW3Pos"
% isW3Pos     = abs(Y0_start)<=100;

%% Initialize the variables
x0_pos      = zeros(Npeaks,1);
y0_pos      = zeros(Npeaks,1);
Sx_pos      = zeros(Npeaks,1);
Sy_pos      = zeros(Npeaks,1);
GSBamp_pos  = zeros(Ndelays,Npeaks);
ESAamp_pos  = zeros(Ndelays,Npeaks);
C_pos       = zeros(Ndelays,Npeaks);
DiagPkXID   = zeros(Npeaks,1);
DiagPkYID   = zeros(Npeaks,1);
pump_pos    = zeros(Npeaks,1); 
probe_pos   = zeros(Npeaks,2); % Probe positions are stored in the form (GSB,ESA)
probe_idx   = zeros(Npeaks,2);

%% Define the fit parameter positions (indexes) of the globally-fit parameters
tot_idx     = 0;
for m=1:Npeaks
    if isDiagonal(m)
        % Define the parameter positions: all diagonal peak parameters are fitted
        x0_pos(m)   = tot_idx+1; tot_idx = tot_idx+1;
        y0_pos(m)   = tot_idx+1; tot_idx = tot_idx+1;
        Sx_pos(m)   = tot_idx+1; tot_idx = tot_idx+1;
        if equal_SxSy == 1
            Sy_pos(m)   = Sx_pos(m);
        else
            Sy_pos(m)   = tot_idx+1; tot_idx = tot_idx+1;
        end
        % Get the starting peak positions (in cm-1) to find later the starting amplitudes
        pump_pos(m)     = X0_start(m);              % W1 position
        probe_pos(m,1)  = X0_start(m);              % W3 position (GSB)
        probe_pos(m,2)  = X0_start(m)-Y0_start(m);  % W3 position (ESA)
    else % If it's not a diagonal peak, X0 and Y0 are given by the corresponding diagonal peaks.
        %%% THIS NEEDS TO BE FIXED PROPERLY FOR INTRAMOLECULAR CROSS PEAKS. PENDING
        % Find the diagonal peak with the closest position
        DiagX               = X0_start;
        DiagX(~isDiagonal)  = NaN;
        DiagPkXID(m)        = findClosestId2Val(DiagX,X0_start(m));
        x0_pos(m)           = x0_pos(DiagPkXID(m));
        % The Y position can be either linked to the diagonal peak or another fit parameter
        if Y0_start(m) < 0
            DiagPkYID(m)        = findClosestId2Val(DiagX,abs(Y0_start(m)));
            y0_pos(m)           = y0_pos(DiagPkYID(m));
            % Get the starting peak positions (in cm-1) to find later the starting amplitudes
            pump_pos(m)     = X0_start(DiagPkXID(m));                           % W1 position
            probe_pos(m,1)  = X0_start(DiagPkYID(m));                           % W3 position (GSB)
            probe_pos(m,2)  = X0_start(DiagPkYID(m))-Y0_start(DiagPkYID(m));    % W3 position (ESA)
        else 
            y0_pos(m) = tot_idx+1; tot_idx = tot_idx+1;
            % Get the starting peak positions (in cm-1) to find later the starting amplitudes
            pump_pos(m)     = X0_start(DiagPkXID(m));                           % W1 position
            probe_pos(m,1)  = X0_start(DiagPkYID(m));                           % W3 position (GSB)
            probe_pos(m,2)  = X0_start(DiagPkYID(m))-Y0_start(DiagPkYID(m));    % W3 position (ESA)
        end
        Sx_pos(m)           = Sx_pos(DiagPkXID(m));
        Sy_pos(m)           = Sy_pos(DiagPkYID(m));

    end
    
end

%% Define the fit parameter positions (indexes) of the amplitudes and correlation coefficients
tot_idx     = max([x0_pos;y0_pos;Sx_pos;Sy_pos]);
peakfnc     = cell(Npeaks,1);

for m=1:Npeaks
    if isDiagonal(m)
        %%% If it's a diagonal peak, 
        for n=1:Ndelays
            GSBamp_pos(n,m) = tot_idx+1; tot_idx = tot_idx+1;
            ESAamp_pos(n,m) = tot_idx+1; tot_idx = tot_idx+1;
            C_pos(n,m)      = tot_idx+1; tot_idx = tot_idx+1;
        end
    else
        %%% If it's not a diagonal peak, X0 and Y0 are given by the corresponding diagonal peaks.
        % Find the diagonal peak with the closest position
        DiagX               = X0_start;
        DiagX(~isDiagonal)  = NaN;
        DiagPkXID(m)        = findClosestId2Val(DiagX,X0_start(m));
        DiagPkYID(m)        = findClosestId2Val(DiagX,abs(Y0_start(m)));
        
        % The positions along W1, Sx and Sy are given by the diagonal peaks
        x0_pos(m) = x0_pos(DiagPkXID(m));
        Sx_pos(m) = Sx_pos(DiagPkXID(m));
        Sy_pos(m) = Sy_pos(DiagPkYID(m));
        
        % If the starting anharmonicity is <0, then link to that of the diagonal peak,
        % otherwise it's a fit parameter
        if Y0_start(m) < 0
            y0_pos(m)   = y0_pos(DiagPkYID(m));
        else
            y0_pos(m)   = tot_idx+1; tot_idx = tot_idx+1;
        end
        
        % The amplitudes are independent for each peak
        for n=1:Ndelays
            GSBamp_pos(n,m) = tot_idx+1; tot_idx = tot_idx+1;
            ESAamp_pos(n,m) = tot_idx+1; tot_idx = tot_idx+1;
        end
        
        % The correlation is either identically zero or another fit parameter
        if ~contains(C_start{m},'0h','IgnoreCase',true)
            for n=1:Ndelays
                C_pos(n,m)  = tot_idx+1; tot_idx = tot_idx+1;
            end
        end
    end
    
    % Build the peak function strings
    x0str   = num2str(x0_pos(m));
    y0str   = num2str(y0_pos(m));
    sxstr   = num2str(Sx_pos(m));
    systr   = num2str(Sy_pos(m));
    peakID  = num2str(m);
    peakGSB = num2str(2.*m-1);
    peakESA = num2str(2.*m);
    switch isDiagonal(m)
        case 1            
            peakfnc{m} = ['G2Dc(X,Y,P(' x0str '),P(' x0str '),P(' sxstr '),P(' systr '),C{' peakID '},A{' peakGSB '}) + G2Dc(X,Y,P(' x0str '),P(' x0str ')-P(' y0str '),P(' sxstr '),P(' systr '),C{' peakID '},A{' peakESA '})'];
        case 0
            % Check if the correlation coefficient of the cross peaks is a fit parameter or if it's set to zero
            if sum(C_pos(:,m),1) == 0
                Cstring = '0';
            else
                Cstring = ['C{' peakID '}'];
            end
            % The position along W1 is given by the diagonal peak
            XpeakGSB_w1 = ['P(' num2str(x0_pos(DiagPkXID(m))) ')'];
            
            % Check if the anharmonicity is the same as the diagonal peak or if it's a fit parameter
            XpeakGSB_w3_str = ['P(' num2str(x0_pos(DiagPkYID(m))) ')']; % GSB of the cross peak between A and B is at (x_A,y_B)
            XpeakESA_w3_str = [XpeakGSB_w3_str '-P(' num2str(y0_pos(DiagPkYID(m))) ')'];
            
            % Write the function string
            peakfnc{m} = ['G2Dc(X,Y,' XpeakGSB_w1 ','  XpeakGSB_w3_str ',P(' sxstr '),P(' systr '),' Cstring ',A{' peakGSB '}) + G2Dc(X,Y,' XpeakGSB_w1 ',' XpeakESA_w3_str ',P(' sxstr '),P(' systr '),' Cstring ',A{' peakESA '})'];
    end
end


%% BUILD THE FIT FUNCTION, STARTING PARAMETERS AND BOUNDS
% Build the fit function string
PeaksFunction     = str2func(['@(P,X,Y,C,A)' strjoin(peakfnc,' + ')]);

% Build the vector of starting parameters
Start_param     = zeros(tot_idx,1); 
UB              = zeros(tot_idx,1); % Upper bounds
LB              = zeros(tot_idx,1); % Lower bounds

%%% Diagonal peaks
% Time-independent parameters
Start_param(x0_pos(isDiagonal)) = X0_start(isDiagonal);
Start_param(y0_pos(isDiagonal)) = Y0_start(isDiagonal);
Start_param(Sx_pos(isDiagonal)) = SX_start(isDiagonal);
Start_param(Sy_pos(isDiagonal)) = SY_start(isDiagonal);

UB(x0_pos(isDiagonal)) = X0_start(isDiagonal)+10;
UB(y0_pos(isDiagonal)) = 50;
UB(Sx_pos(isDiagonal)) = 50;
UB(Sy_pos(isDiagonal)) = 50;

LB(x0_pos(isDiagonal)) = X0_start(isDiagonal)-10;
LB(y0_pos(isDiagonal)) = 5;
LB(Sx_pos(isDiagonal)) = 5;
LB(Sy_pos(isDiagonal)) = 5;

%%% Cross peaks
%%%%% PENDING: Explicit definition for intramolecular cross peaks

%%% Time-dependent parameters
% Spectral diffusion
Start_param(C_pos(C_pos~=0)) = 0.4;     % Starting C value of 0.5 for all points. Needs improvement.
UB(C_pos(C_pos~=0))     = 0.99;         % C can only be in the range [-1 1]. Considering 0 to 1 only
LB(C_pos(C_pos~=0))     = 0;            % C can only be in the range [-1 1]. Considering 0 to 1 only


% Amplitudes
% The data is normalized, so the min/max amplitudes should be around +- 5
    % Get the positions of each peak in terms of the W1 and W3 axes of the dataset
    pump_idx        = findClosestId2Val(PumpAxis,pump_pos)';
    probe_idx(:,1)  = findClosestId2Val(ProbeAxis,probe_pos(:,1));
    probe_idx(:,2)  = findClosestId2Val(ProbeAxis,probe_pos(:,2));

    % Evaluate the amplitudes of the input data at the position of each peak and use them as initial
    % parameters for the fit
    for m=1:Npeaks
        Start_param(GSBamp_pos(:,m))  = squeeze(ZData(pump_idx(m),probe_idx(m,1),:)); % GSB
        Start_param(ESAamp_pos(:,m))  = squeeze(ZData(pump_idx(m),probe_idx(m,2),:))'; % ESA
    end

UB(GSBamp_pos(:))       = 0;            % Bleaches can only be negative
UB(ESAamp_pos(:))       = +5;           % ESA can only be positive

LB(GSBamp_pos(:))       = -5;           % Bleaches can only be negative
LB(ESAamp_pos(:))       = 0;            % ESA can only be positive


%% SAVE the positions of the parameters to the ParamPos structure
ParamPos.Npeaks     = Npeaks;
ParamPos.Ndelays    = Ndelays;
ParamPos.x0_pos     = x0_pos;
ParamPos.y0_pos     = y0_pos;
ParamPos.Sx_pos     = Sx_pos;
ParamPos.Sy_pos     = Sy_pos;
ParamPos.GSBamp_pos = GSBamp_pos;
ParamPos.ESAamp_pos = ESAamp_pos;
ParamPos.isDiagonal = isDiagonal;
ParamPos.C_pos      = C_pos;