function stop = optimplotx_txt(x,optimValues,state,varargin)
% OPTIMPLOTX Plot current point at each iteration.
%
%   STOP = OPTIMPLOTX(X,OPTIMVALUES,STATE) plots the current point, X, as a
%   bar plot of its elements at the current iteration.
%
%   Example:
%   Create an options structure that will use OPTIMPLOTX
%   as the plot function
%       options = optimset('PlotFcns',@optimplotx);
%
%   Pass the options into an optimization problem to view the plot
%       fminbnd(@sin,3,10,options)

%   Copyright 2006-2010 The MathWorks, Inc.
% MODIFIED BY: Ricardo J. Fernández-Terán

stop = false;
switch state
    case 'iter'
        % Reshape if x is a matrix
        x = x(:);
        xLength = length(x);
        xlabelText = getString(message('MATLAB:optimfun:funfun:optimplots:LabelNumberOfVariables',sprintf('%g',xLength)));

        % Display up to the first 100 values
        if length(x) > 100
            x = x(1:100);
            xlabelText = {xlabelText,getString(message('MATLAB:optimfun:funfun:optimplots:LabelShowingOnlyFirst100Variables'))};
        end
        
        

        if optimValues.iteration == 0
            % The 'iter' case is  called during the zeroth iteration,
            % but it now has values that were empty during the 'init' case
            plotx = bar(x);
            title(getString(message('MATLAB:optimfun:funfun:optimplots:TitleCurrentPoint')),'interp','none');
            ylabel(getString(message('MATLAB:optimfun:funfun:optimplots:LabelCurrentPoint')),'interp','none');
            xlabel(xlabelText,'interp','none');
            set(plotx,'edgecolor','none')
            set(gca,'xlim',[0,1 + xLength])
            set(plotx,'Tag','optimplotx');
        else
            plotx = findobj(get(gca,'Children'),'Tag','optimplotx');
            set(plotx,'Ydata',x);
            delete(findobj(get(gca,'Children'),'Type','Text'));
        end
        xtips1 = plotx(1).XEndPoints;
        ytips1 = plotx(1).YEndPoints;
        labels1 = num2str(plotx(1).YData','%.3f');
        text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom');
end

