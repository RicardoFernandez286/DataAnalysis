function L2Ddata = L2Dc(X,Y,x0,y0,sx,sy,c,varargin)
    % Description: Evaluates a two-dimensional Gaussian function with correlation
    % over the XY mesh given by X and Y, with center coordinates (x0,y0) and
    % standard deviations sx and sy. The parameter c is a correlation
    % coefficient, related to the tilt of the curve along the diagonal,
    % where 0 is a normal round Gaussian and ~1 gives a diagonal Gaussian.
    %
    % Inputs:
    %   X,Y     = Mesh over which the function will be evaluated. 
    %             To create a mesh from known linear intervals, use the "meshgrid" function.
    %   x0,y0   = Centroids of the Lorentzian function.
    %   sx,sy   = Sigma values along the X and Y dimensions, respectively.
    %   c       = Correlation coefficient ranging from -1 to 1 (typically positive).
    %
    % Ricardo Fernández-Terán / 2019.03.30 / v1.0b
    
    % Check if the amplitude is given or not (set it to 1 if not)
    if isempty(varargin)
        A = 1;
    else
        A = varargin{1};
    end
    
    % Calculate the output
    calc = (A./(1-c.^2).^(1/2)).*(1./(1 + ((X-x0)./sx).^2 + ((Y-y0)./sy).^2 - (2*c.*(X-x0).*(Y-y0))./(sx.*sy)).^(3/2));
    
    % If operating in 4D, reduce the output to take only the elements with indexes XYZZ
    % (to have A(t) and C(t) at the same t)
    if length(A) > 1 || length(c) > 1
        if size(calc,3) == size(calc,4)
            Size4D  = size(calc);
            L2Ddata     = zeros(Size4D(1:3));
            for q=1:Size4D(3)
                L2Ddata(:,:,q) = calc(:,:,q,q);
            end
        else
            error('Incosistent A(t) and C(t) parameter sizes')
        end
    else
        L2Ddata=calc;
    end
end