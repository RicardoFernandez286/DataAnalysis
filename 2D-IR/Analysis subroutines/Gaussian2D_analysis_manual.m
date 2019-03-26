function fitted_param = Gaussian2D_analysis_manual
clear all

Npix        = 32;
Nw1         = 70;
Ncontours   = 40;
PlotData    = 1;
PlotSimu    = 0;
PlotDiff    = 1;
DoFit       = 1;

if PlotSimu == 1
    %% Plot a 2D Gaussian function over a grid
    % Define a grid
    w1       = linspace(1920,2120,Nw1);
    w3       = linspace(1920,2120,Npix);

    [W1,W3] = meshgrid(w1,w3);

    % Calculate the peaks
    AllPeaks = zeros(Npix,Nw1,8);

    AllPeaks(:,:,1) = Gaussian2D(W1,W3,1985,1985,20,12,0.0,-1);
    AllPeaks(:,:,2) = Gaussian2D(W1,W3,1985,1960,20,12,0.0,+0.8);
    AllPeaks(:,:,3) = Gaussian2D(W1,W3,2062,2062,20,12,0.0,-1);
    AllPeaks(:,:,4) = Gaussian2D(W1,W3,2062,2037,20,12,0.0,+0.8);

    AllPeaks(:,:,5) = 0.5*Gaussian2D(W1,W3,1985,2062,20,12,0.0,-1);
    AllPeaks(:,:,6) = 0.5*Gaussian2D(W1,W3,1985,2037,20,12,0.0,+0.8);
    AllPeaks(:,:,7) = 0.5*Gaussian2D(W1,W3,2062,1985,20,12,0.0,-1);
    AllPeaks(:,:,8) = 0.5*Gaussian2D(W1,W3,2062,1960,20,12,0.0,+0.8);

    AllPeaks = sum(AllPeaks,3);
    toc
    contourf(w1,w3,AllPeaks,Ncontours,'LineColor','flat','LineStyle','-','LineWidth',0.5)
    colorbar
    % Change the colormap to DkRed/White/DkBlue
    cmap            = darkb2r(-1,1,50,2);
    colormap(cmap)

    diagline        = refline(1,0);
    diagline.Color  = [0 0 0];
end

%% Read the data
data    = dlmread('SolC_600fs.dat');
inputW1 = data(1,2:end)';
inputW3 = data(2:end,1);
Zdata   = data(2:end,2:end);
Zdata   = Zdata./max(Zdata(:));

Omega = {inputW1;inputW3};

%% DO the fit
if DoFit == 0
    return
end
%         w01  w02 D1 D2 sx sy  C  ab1 ab2 ae1 ae2
% testP = [1985 2062 25 25 20 12 0.8 1   1   0.8 0.8  ];
% LB    = [1960 2050 5  5  5  5  0   0   0   0   0    ];
% UB    = [2000 2075 50 50 30 30 1   5   5   5   5    ];

testP = [1940 2030 20 20 5 5   0.5  1   1   1   0    0   0   0   0   ];
LB    = [1900 2020 5  5  5  5  0    0   0   0   0    0   0   0   0   ];
UB    = [2000 2100 50 50 30 30 1  	9   9   9   9    9   9   9   9  ];

options = optimoptions('lsqcurvefit',...
            'MaxFunctionEvaluations',5000,...
            'MaxIterations',5000,...
            'Algorithm','trust-region-reflective',...  %'levenberg-marquardt' 'trust-region-reflective'
            'OptimalityTolerance',1e-15,...
            'FunctionTolerance',1e-15,...
            'StepTolerance',1e-15);

tic        
[fitted_param,resnorm,residuals,exitflag,output,lambda,jacobian_fit] = lsqcurvefit(@FitFunction,testP,Omega,Zdata,LB,UB,options);
toc

fitted_param = fitted_param';

%% Plot the fit result
fh = figure(1);
clf(fh)
ax1 = axes('parent',fh);
contourf(ax1,inputW1,inputW3,FitFunction(fitted_param,Omega),Ncontours,'LineColor','flat','LineStyle','-','LineWidth',0.5)
colorbar
% Change the colormap to DkRed/White/DkBlue
cmap            = darkb2r(-1,1,Ncontours,2);
colormap(cmap)
caxis([-1 1])
diagline        = refline(1,0);
diagline.Color  = [0 0 0];

%% Plot the original data
if PlotData == 1
    hold(ax1,'on')
    contour(ax1,inputW1,inputW3,Zdata,Ncontours,'LineColor',[0.8 0.8 0.8],'LineStyle','-','LineWidth',1)
    colorbar
    diagline        = refline(1,0);
    diagline.Color  = [0 0 0];
    caxis([-1 1])
end

%% Plot the difference
if PlotDiff == 1
    fh = figure(2);
    clf(fh)
    ax2 = axes('parent',fh);
    contourf(ax2,inputW1,inputW3,Zdata-FitFunction(fitted_param,Omega),Ncontours,'LineColor','flat','LineStyle','-','LineWidth',0.5)
    colorbar
    % Change the colormap to DkRed/White/DkBlue
    cmap            = darkb2r(-0.1,0.1,Ncontours,2);
    colormap(cmap)
    caxis(ax2,[-0.1 0.1])
    diagline        = refline(1,0);
    diagline.Color  = [0 0 0];
end

function output = FitFunction(P,XY)
    % Fit function with parameters P and XY
    % P is a vector of fit parameters containing the following things:
        % P(1) = bleach1
        % P(2) = bleach2
        % P(3) = Anharm1
        % P(4) = Anharm2
        % P(5) = sigma_x
        % P(6) = sigma_y
        % P(7) = C
        % P(8) = A(bleach1)
        % P(9) = A(bleach2)
        % P(10)= A(ESA1)
        % P(11)= A(ESA2)
        % P(12)= A(upperXP_b)
        % P(13)= A(lowerXP_b)
        % P(14)= A(upperXP_e)
        % P(15)= A(lowerXP_e)
    % XY is a 2-element cell containing the elements of the W1 and W3 axes, in that order.
    
    % Create a mesh grid to evaluate the fit function
    [W1_fit,W3_fit] = meshgrid(XY{1},XY{2});
    
    % Build the peaks
    Peaks = zeros(length(XY{2}),length(XY{1}),8);
    for i=1:2
        Peaks(:,:,2*i-1)    = -Gaussian2D(W1_fit,W3_fit,P(i),P(i),P(5),P(6),P(7),P(7+i));       % Bleach
        Peaks(:,:,2*i)      = Gaussian2D(W1_fit,W3_fit,P(i),P(i)-P(2+i),P(5),P(6),P(7),P(9+i)); % ESA
    end
    Peaks(:,:,5)    = -Gaussian2D(W1_fit,W3_fit,P(1),P(2),P(5),P(6),0,P(12));       % Bleach Xpeak1
    Peaks(:,:,6)    = Gaussian2D(W1_fit,W3_fit,P(1),P(2)-P(4),P(5),P(6),0,P(14)); % ESA Xpeak1
    Peaks(:,:,7)    = -Gaussian2D(W1_fit,W3_fit,P(2),P(1),P(5),P(6),0,P(13));       % Bleach Xpeak2
    Peaks(:,:,8)    = Gaussian2D(W1_fit,W3_fit,P(2),P(1)-P(3),P(5),P(6),0,P(15)); % ESA Xpeak2        
    
    output = sum(Peaks,3);
    
end

end