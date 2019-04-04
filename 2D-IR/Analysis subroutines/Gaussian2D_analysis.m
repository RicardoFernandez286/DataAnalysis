function handles = Gaussian2D_analysis(handles)
% Description: This function performs a 2D Gaussian fitting analysis of a
% series of 2D-IR data as a function of the waiting time
% 
% Usage: handles = Gaussian2D_analysis(handles)
% Inputs: handles structure with the following fields:
%     ProbeAxis
%     PumpAxis
%     t2delays
%     PROC_2D_DATA
%     InteractiveModeTick.Value (a tick box to define whether the user is asked for graphical or textual input)
% Outputs:
%     fitPar and firErr: two data structures containing the fit results and their standard errors.
%
% Ricardo Fernández-Terán / 04.04.2019 / v2.2a

%% READ from handles
% Get information from the GUI
ProbeAxis           = handles.ProbeAxis;
PumpAxis            = handles.PumpAxis{1,1};
t2delays            = handles.t2delays;
PROC_2D_DATA        = handles.PROC_2D_DATA;

% Hardcoded settings
% Minimum t2 time to start the fit (will remove all t2 delays before)
t2_startFit         = 0.1; % in ps

% Initialize variables (ensures only the current analysis is saved)


%% Cut the data
cut_data = questdlg('Cut the 2D dataset?','Cut 2D dataset','Yes, define region','Use probe axis','No','Use probe axis');

switch cut_data
    case 'Yes, define region'
        switch handles.InteractiveModeTick.Value
            case 1
                % Select points interactively, to finish hit RETURN
                [RegionPositions,~] = SelectRegions(Axes);
                % Separate the values into pump and probe, according to how the data is plotted
                % i.e. convert from [x_min x_max y_min y_max] to [pump_min pump_max probe_min probe_max]
                % Get only ONE region
                switch plot_pumpdirection
                    case 'Horizontal'
                        pump_range      = RegionPositions(1,1:2); 
                        probe_range     = RegionPositions(1,3:4);
                    case 'Vertical'
                        pump_range      = RegionPositions(1,3:4); 
                        probe_range     = RegionPositions(1,1:2);
                end
            case 0
                RegionPositions_cell    = inputdlg({'Enter the pump region (Format: minWL maxWL): ','Enter the probe regions (Format: minWL maxWL):'},'Define regions to integrate', [1 80]);
                pump_range              = str2num(RegionPositions_cell{1});
                probe_range             = str2num(RegionPositions_cell{2});
                % Get only 1 region    
                pump_range              = pump_range(1);
                probe_range             = probe_range(1);
                % Check that everything is in order
                L = size(pump_range,1);
                M = size(probe_range,1);
                if L ~= M
                    warndlg('Number of pump and probe ranges is different!')
                    return
                end
        end
    case 'Use probe axis'
        probe_range = [min(ProbeAxis) max(ProbeAxis)];
        pump_range  = probe_range;
    case 'No'
        probe_range = [min(ProbeAxis) max(ProbeAxis)];
        pump_range  = [min(PumpAxis) max(PumpAxis)];
end

% Get the indices
pump_idxrange   = sort(findClosestId2Val(PumpAxis,pump_range));
probe_idxrange  = sort(findClosestId2Val(ProbeAxis,probe_range));

% Consider only delays after a given t2 delay
PROC_2D_DATA    = PROC_2D_DATA(t2delays>t2_startFit,:);
t2delays        = t2delays(t2delays>t2_startFit);
Ndelays         = length(t2delays);
Znormfactor     = zeros(Ndelays,1);
ZData           = zeros(pump_idxrange(2)-pump_idxrange(1)+1,probe_idxrange(2)-probe_idxrange(1)+1,Ndelays);

% Cut the data
for m=1:Ndelays
    data                = PROC_2D_DATA{m,1};
    PROC_2D_DATA{m,1}   = data(pump_idxrange(1):pump_idxrange(2),probe_idxrange(1):probe_idxrange(2));
    Znormfactor(m)      = max(abs(PROC_2D_DATA{m,1}),[],'all');
    ZData(:,:,m)        = PROC_2D_DATA{m,1}./Znormfactor(m);
end

% Cut the pump and probe axes
PumpAxis            = PumpAxis(pump_idxrange(1):pump_idxrange(2));
ProbeAxis           = ProbeAxis(probe_idxrange(1):probe_idxrange(2));
Omega               = {PumpAxis;ProbeAxis};

%% Get the user's starting parameters - v1.0, without preview
% Get the parameters
fitparameters = Gaussian2D_fitparam;

% If the user cancelled, return
if isempty(fitparameters)
    return
end

% Parse the input into a fit function, calculate starting parameters and limits for the fit
[PeaksFunction,Start_param,UB,LB,ParamPos] = parse_2DGC_input(fitparameters,Ndelays,Omega,ZData);

%% DO the fit
options = optimoptions('lsqcurvefit',...
            'MaxFunctionEvaluations',6000,...
            'MaxIterations',50,...
            'Algorithm','trust-region-reflective',...  %'levenberg-marquardt' 'trust-region-reflective'
            'OptimalityTolerance',5e-5,...
            'FunctionTolerance',5e-5,...
            'StepTolerance',5e-5,...
            'UseParallel',true,...
            'TypicalX',Start_param,...
            'SubproblemAlgorithm','factorization',...
            'PlotFcn','optimplotresnorm');
        
% options = optimoptions('lsqcurvefit',...
%             'MaxFunctionEvaluations',5000,...
%             'MaxIterations',10,...
%             'Algorithm','levenberg-marquardt',...  %'levenberg-marquardt' 'trust-region-reflective'
%             'OptimalityTolerance',1e-5,...
%             'FunctionTolerance',1e-5,...
%             'StepTolerance',1e-4,...
%             'PlotFcn','optimplotresnorm',...
%             'UseParallel',true,...
%             'TypicalX',Start_param);
% %             'ScaleProblem','jacobian',...

% Build the input structure with the necessary ingredients
input.Omega         = Omega;
input.ParamPos      = ParamPos;
input.PkFunction    = PeaksFunction;

% Perform the actual fit
t_start = tic;
% [fitted_param,resnorm,residuals,exitflag,output_st,lambda,jacobian_fit] = lsqcurvefit(@FitFunction,Start_param,input,ZData,[],[],options);
[fitted_param,SSR,residuals,exitflag,output_st,~,jacobian_fit] = lsqcurvefit(@FitFunction,Start_param,input,ZData,LB,UB,options);
t_fit   = toc(t_start);

if exitflag >= 0
    msgbox({['Fit completed in ' num2str(t_fit,'%.4g') ' seconds.'];[];['SSR = ' num2str(SSR,'%.5g')];['Iterations = ' num2str(output_st.iterations)];['Func. evals. = ' num2str(output_st.funcCount)];['Step size = ' num2str(output_st.stepsize)]},'Fit completed!','help');
else
    msgbox('Error in the fit! Please check all parameters','Error during fit','error');
end

%% Evaluate the solution
FitResults = bsxfun(@times,FitFunction(fitted_param,input),reshape(Znormfactor,1,1,[]));

% Get the solution parameters
fitPar.X0           = fitted_param(ParamPos.x0_pos);
fitPar.Anharm       = fitted_param(ParamPos.y0_pos);
fitPar.Y0           = fitted_param(ParamPos.x0_pos)-fitted_param(ParamPos.y0_pos);
fitPar.Sx           = fitted_param(ParamPos.Sx_pos);
fitPar.Sy           = fitted_param(ParamPos.Sy_pos);
fitPar.C            = fitted_param(ParamPos.C_pos(:,ParamPos.isDiagonal));
fitPar.Amps(:,:,1)  = fitted_param(ParamPos.GSBamp_pos).*Znormfactor;
fitPar.Amps(:,:,2)  = fitted_param(ParamPos.ESAamp_pos).*Znormfactor;

% Calculate the errors of the fitted parameters
fitPar_CI = nlparci(fitted_param,residuals,'jacobian',jacobian_fit);
fitted_err= (fitPar_CI(:,2)-fitPar_CI(:,1))./2;

fitErr.X0           = fitted_err(ParamPos.x0_pos);
fitErr.Anharm       = fitted_err(ParamPos.y0_pos);
fitErr.Y0           = fitted_err(ParamPos.x0_pos)+fitted_err(ParamPos.y0_pos);
fitErr.Sx           = fitted_err(ParamPos.Sx_pos);
fitErr.Sy           = fitted_err(ParamPos.Sy_pos);
fitErr.C            = fitted_err(ParamPos.C_pos(:,ParamPos.isDiagonal));
fitErr.Amps(:,:,1)  = fitted_err(ParamPos.GSBamp_pos).*Znormfactor;
fitErr.Amps(:,:,2)  = fitted_err(ParamPos.ESAamp_pos).*Znormfactor;

% Calculate the amplitude to volume conversion factor (including C)
corrterm                        = zeros(Ndelays,ParamPos.Npeaks);
corrterm(:,ParamPos.isDiagonal) = fitted_param(ParamPos.C_pos(:,ParamPos.isDiagonal));
corrterm                        = sqrt(1-corrterm.^2);

fitPar.Vols(:,:,1)  = fitPar.Amps(:,:,1).*corrterm;
fitPar.Vols(:,:,2)  = fitPar.Amps(:,:,2).*corrterm;
fitErr.Vols(:,:,1)  = fitErr.Amps(:,:,1).*corrterm;
fitErr.Vols(:,:,2)  = fitErr.Amps(:,:,2).*corrterm;

% Normalize the data in the "usual" way
NormVols(:,[1 3],1) = -fitPar.Vols(:,[1 3],1)./max(abs(fitPar.Vols(:,1,1)));
NormVols(:,[1 3],2) = fitPar.Vols(:,[1 3],2)./max(abs(fitPar.Vols(:,1,2)));
NormVols(:,[2 4],1) = -fitPar.Vols(:,[2 4],1)./max(abs(fitPar.Vols(:,2,1)));
NormVols(:,[2 4],2) = fitPar.Vols(:,[2 4],2)./max(abs(fitPar.Vols(:,2,2)));

NormErr(:,[1 3],1)  = fitErr.Vols(:,[1 3],1)./max(abs(fitPar.Vols(:,1,1)));
NormErr(:,[1 3],2)  = fitErr.Vols(:,[1 3],2)./max(abs(fitPar.Vols(:,1,2)));
NormErr(:,[2 4],1)  = fitErr.Vols(:,[2 4],1)./max(abs(fitPar.Vols(:,2,1)));
NormErr(:,[2 4],2)  = fitErr.Vols(:,[2 4],2)./max(abs(fitPar.Vols(:,2,2)));

NormVols(:,[3 4],:) = 10*NormVols(:,[3 4],:);
NormErr(:,[3 4],:)  = 10*NormErr(:,[3 4],:);

% % dlmwrite('SolC_fittedC.dat',[t2delays fitPar.C fitErr.C]);
% dlmwrite('SolC_fittedAmps_GSB.dat',[t2delays NormAmps(:,:,1) NormErr(:,:,1)]);
% dlmwrite('SolC_fittedAmps_ESA.dat',[t2delays NormAmps(:,:,2) NormErr(:,:,2)]);
% 
% %% Plot the results :)
% 
% Plot the fit results at a given delay
% N=19;
% figure(3)
% contourf(Omega{1},Omega{2},FitResults(:,:,N)',40,'LineColor','flat');
% hold on
% contour(Omega{1},Omega{2},ZData(:,:,N)',40,'LineColor',0.8*[1 1 1]);
% hold off
% % contourf(Omega{1},Omega{2},(ZData(:,:,N)-FitResults(:,:,N))',30,'LineColor','flat'); 
% colorbar; 
% colormap(darkb2r(-1,1,80,2)); 
% diagline = refline(1,0); diagline.Color  = [0 0 0];
% 
% pepita=1;
% 
% figure; plot(t2delays,-fitPar.Amps(:,[1 3],1))

%% Save the fit results to a file
filename    = char([handles.rootdir filesep char(handles.datafilename) '_FIT_RESULTS.mat']);
save(filename,'fitPar','fitErr','SSR','output_st','t2delays','input','FitResults','NormVols','NormErr');

%% Save the fit to handles
handles.FitResults  = FitResults;
handles.t2_startFit = t2_startFit;
handles.FitInput    = input;

%% FUNCTION DEFINITIONS

% Build the fit function
function output3D = FitFunction(P,Input)
    W1W3_Axes   = Input.Omega;
    Param_Pos   = Input.ParamPos;
    PkFunction  = Input.PkFunction;
    C   = cell(Param_Pos.Npeaks,1);
    A   = cell(2*Param_Pos.Npeaks,1);
    for k=1:Param_Pos.Npeaks
        if sum(Param_Pos.C_pos(:,k)) ~= 0
            [~,~,~   ,A{2*k-1}] = ndgrid(W1W3_Axes{1},W1W3_Axes{2},P(Param_Pos.C_pos(:,k)),P(Param_Pos.GSBamp_pos(:,k))); % GSB
            [X,Y,C{k},A{2*k}]   = ndgrid(W1W3_Axes{1},W1W3_Axes{2},P(Param_Pos.C_pos(:,k)),P(Param_Pos.ESAamp_pos(:,k))); % ESA
        else
            C{k} = zeros(Param_Pos.Ndelays,1);
            [~,~,~   ,A{2*k-1}] = ndgrid(W1W3_Axes{1},W1W3_Axes{2},zeros(Param_Pos.Ndelays,1),P(Param_Pos.GSBamp_pos(:,k))); % GSB
            [X,Y,C{k},A{2*k}]   = ndgrid(W1W3_Axes{1},W1W3_Axes{2},zeros(Param_Pos.Ndelays,1),P(Param_Pos.ESAamp_pos(:,k))); % ESA
        end
    end
    % Evaluate the fit function with the given parameters
    output3D = PkFunction(P,X,Y,C,A);
end

end