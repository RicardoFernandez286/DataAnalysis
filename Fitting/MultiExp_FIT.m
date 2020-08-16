t       = t2delays;
ZData   = ESA_down;


Klb     = [60,  150,   70*8   ];
K       = [70,  70*8,  70*64  ];
Kub     = [85,  70*16, 70*64*2];

Ntraces = size(ZData,2);
Ntau    = length(K);


A           = ones(Ntraces,Ntau+1);
Alb         = -A*Inf;
Aub         = A*0;
Aub(:,1)    = A(:,1)*70;
Alb(:,1)    = A(:,1)*40;

x0  = [K(:); A(:)];

LB  = [Klb(:); Alb(:)];
UB  = [Kub(:); Aub(:)];

options = optimoptions('lsqcurvefit',...
                    'MaxFunctionEvaluations',6500,...
                    'MaxIterations',250,...
                    'Algorithm','trust-region-reflective',...  %'levenberg-marquardt' 'trust-region-reflective'
                    'OptimalityTolerance',1e-10,...
                    'FunctionTolerance',eps,...
                    'StepTolerance',eps,...
                    'UseParallel',true,...
                    'SubproblemAlgorithm','factorization',...
                    'PlotFcn','optimplotresnorm');
                
% Perform the actual fit
    t_start = tic;
    % [fitted_param,resnorm,residuals,exitflag,output_st,lambda,jacobian_fit] = lsqcurvefit(@FitFunction,Start_param,input,ZData,[],[],options);
    [fitted_param,SSR,residuals,exitflag,output_st,~,jacobian_fit] = lsqcurvefit(@(P,t)FitFunc(P,t,Ntraces,Ntau),x0,t,ZData,LB,UB,options);
    t_fit   = toc(t_start);

Zfit = FitFunc(fitted_param,t,Ntraces,Ntau);

%Calculate parameters and errors
fitted_CI = nlparci(fitted_param,residuals,'jacobian',jacobian_fit);
fitted_SD = fitted_CI(:,2) - fitted_CI(:,1);

Kfit = fitted_param(1:Ntau);
Afit = reshape(fitted_param((Ntau+1):end),[Ntraces,Ntau+1]);

Kfit_SD = fitted_SD(1:Ntau);
Afit_SD = reshape(fitted_SD((Ntau+1):end),[Ntraces,Ntau+1]);

% Plot fit results
cmap = othercolor('Mrainbow',Ntraces);
fh = figure(1); clf(fh);
ax = axes('parent',fh);
hold(ax,'on');

for i=1:Ntraces
    plot(ax,t,ZData(:,i),'o','Color',cmap(i,:),'MarkerSize',6,'HandleVisibility','off');
    plot(ax,t,Zfit(:,i),'-','Color',cmap(i,:),'LineWidth',2,'DisplayName',num2str(100-ConcPercent(i)));
end
hold(ax,'off');

ax.FontSize = 16;
axis(ax,'tight');

xlabel(ax,'t_{2} delay (ps)','FontWeight','bold');
ylabel(ax,'Population (%)','FontWeight','bold');

box(ax,'on');
lg = legend(ax,'show','location','best');
lg = legend(ax,'boxoff');
title(lg,'% CNBz');

%%
function Y = FitFunc(x,t,Ntraces,Ntau)
    K = x(1:Ntau);
    A = reshape(x((Ntau+1):end),[Ntraces,Ntau+1]);
    Y = MultiExp(K,A,t);
end

function Y = MultiExp(K,A,t)
   Y = A(:,1).*ones(1,length(t)) + A(:,2:end)*exp(-t'./K);
   Y = Y';
end