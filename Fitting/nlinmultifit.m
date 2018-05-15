function [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit(x_cell, y_cell, mdl_cell, beta0, varargin)
%NLINMULTIFIT Nonlinear least-squares regression of multiple data sets
%
%	A wrapper function for NLINFIT which allows simulatenous fitting of
%	multiple data sets with shared fitting parameters. See example below.
%
%	Unlike different solutions (using fminsearch or similar functions)
%	this approach enables simple estimation of model predictions and
%	their confidence intervals, as well as confidence intervals on the
%	fitted parameters (using the built-in NLPREDCI and NLPARCI functions).
%
%	INPUT:
% 		x_cell,y_cell: Cell arrays containing the x,y vectors of the fitted
%					   data sets.
% 		mdl_cell: Cell array containing model functions for each data set.
% 		beta0: Vector containing initial guess of the fitted parameters.
% 		options: Structure containing control parameters for NLINFIT (see
%				 help file on NLINFIT for more details).
%		Name,Value pairs: additional options passed to NLINFIT (see help
%						  file on NLINFIT for more details). The 'Weights'
%						  input option may be specified as a cell array,
%						  similar to the input of x_cell and y_cell.
%
%	NOTE:
%		Unless a cell array of weights is explicitly given, NLINMULTIFIT
%		will create a weights vector such that all points in all datasets
%		receive equal weight regardless of each dataset length.
%		
%
%	OUTPUT:
%		beta,r,J,Sigma,mse,errorparam,robustw: Direct output from NLINFIT.
%
%	EXAMPLE:
% 		% Generate X vectors for both data sets
% 		x1 = 0:0.1:10;
% 		x2 = 0:1:10;
% 
% 		% Generate Y data with some noise
% 		y1 = cos(2*pi*0.5*x1).*exp(-x1/5) + 0.05*randn(size(x1));
% 		y2 = 0.5 + 2*exp(-(x2/5)) + 0.05*randn(size(x2));
% 
% 		% Define fitting functions and parameters, with identical
%		% exponential decay for both data sets
% 		mdl1 = @(beta,x) cos(2*pi*beta(1)*x).*exp(-x/beta(2));
% 		mdl2 = @(beta,x) beta(4) + beta(3)*exp(-(x/beta(2)));
% 
% 		% Prepare input for NLINMULTIFIT and perform fitting
% 		x_cell = {x1, x2};
% 		y_cell = {y1, y2};
% 		mdl_cell = {mdl1, mdl2};
% 		beta0 = [1, 1, 1, 1];
% 		[beta,r,J,Sigma,mse,errorparam,robustw] = ...
%					nlinmultifit(x_cell, y_cell, mdl_cell, beta0);
% 
% 		% Calculate model predictions and confidence intervals
% 		[ypred1,delta1] = nlpredci(mdl1,x1,beta,r,'covar',Sigma);
% 		[ypred2,delta2] = nlpredci(mdl2,x2,beta,r,'covar',Sigma);
% 
% 		% Calculate parameter confidence intervals
% 		ci = nlparci(beta,r,'Jacobian',J);
% 
% 		% Plot results
% 		figure;
% 		hold all;
% 		box on;
% 		scatter(x1,y1);
% 		scatter(x2,y2);
% 		plot(x1,ypred1,'Color','blue');
% 		plot(x1,ypred1+delta1,'Color','blue','LineStyle',':');
% 		plot(x1,ypred1-delta1,'Color','blue','LineStyle',':');
% 		plot(x2,ypred2,'Color',[0 0.5 0]);
% 		plot(x2,ypred2+delta2,'Color',[0 0.5 0],'LineStyle',':');
% 		plot(x2,ypred2-delta2,'Color',[0 0.5 0],'LineStyle',':');
%
%	AUTHOR:
%		Chen Avinadav
%		mygiga (at) gmail
%
		
	if (~iscell(x_cell) && ~iscell(y_cell) && ~iscell(mdl_cell))
		[beta,r,J,Sigma,mse,errorparam,robustw] = nlinfit(x_cell, y_cell, mdl_cell, beta0, varargin{1:end});
		return;
	end
	
	num_curves = length(x_cell);
	if length(y_cell) ~= num_curves || length(mdl_cell) ~= num_curves
		error('Invalid input to NLINMULTIFIT');
	end
	
	x_vec = [];
	y_vec = [];
	mdl_vec = '@(beta,x) [';
	mdl_ind1 = 1;
	mdl_ind2 = 0;
	
	found_weights = 0;
	for ii = 1:length(varargin)
		if strcmpi(varargin{ii}, 'weights')
			found_weights = 1;
			break;
		end
	end
	if ~found_weights
		wghts_cell = {};
		for ii = 1:num_curves
			wghts_cell{ii} = @(y) ones(size(y))/length(x_cell{ii});
		end
		varargin{end+1} = 'Weights';
		varargin{end+1} = wghts_cell;
	end
	
	wghts_argin_ind = 0;
	wghts_cell = {};
	for ii = 1:length(varargin)
		if strcmpi(varargin{ii}, 'weights')
			wghts_argin_ind = ii+1;
			wghts_cell = varargin{ii+1};
			if isa(wghts_cell{1}, 'function_handle')
				is_wght_func = 1;
				wghts_vec = '@(y) [';
			else
				is_wght_func = 0;
				wghts_vec = [];
			end
			break;
		end
	end
	
	for ii = 1:num_curves
		if length(x_cell{ii}) ~= length(y_cell{ii})
			error('Invalid input to NLINMULTIFIT');
		end
		if size(x_cell{ii},2) == 1
			x_cell{ii} = x_cell{ii}';
		end
		if size(y_cell{ii},2) == 1
			y_cell{ii} = y_cell{ii}';
		end
		x_vec = [x_vec, x_cell{ii}];
		y_vec = [y_vec, y_cell{ii}];
		mdl_ind2 = mdl_ind2 + length(x_cell{ii});
		mdl_vec = [mdl_vec, sprintf('mdl_cell{%d}(beta,x(%d:%d)), ', ii, mdl_ind1, mdl_ind2)];
		if ~isempty(wghts_cell)
			if is_wght_func
				wghts_vec = [wghts_vec, sprintf('wghts_cell{%d}(y(%d:%d)), ', ii, mdl_ind1, mdl_ind2)];
			else
				if size(wghts_cell{ii},2) == 1
					wghts_cell{ii} = wghts_cell{ii}';
				end
				wghts_vec = [wghts_vec, wghts_cell{ii}];
			end
		end
		mdl_ind1 = mdl_ind1 + length(x_cell{ii});
	end
	mdl_vec = [mdl_vec(1:end-2), '];'];
	mdl_vec = eval(mdl_vec);
	if ~isempty(wghts_cell)
		if is_wght_func
			wghts_vec = [wghts_vec(1:end-2), '];'];
			wghts_vec = eval(wghts_vec);
		end
		varargin{wghts_argin_ind} = wghts_vec;		
	end
		
	[beta,r,J,Sigma,mse,errorparam,robustw] = nlinfit(x_vec, y_vec, mdl_vec, beta0, varargin{1:end});
end
