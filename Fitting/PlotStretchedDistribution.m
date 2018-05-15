function PlotStretchedDistribution(tauKWW,beta,Nterms,timescale,tmin,tmax,plotscale,ShowTauAvg,indexing)
% Syntax: PlotStretchedDistribution(tauKWW,beta,Nterms,timescale,tmin,tmax,plotscale)
% Plots the Kohlrausch-Williams-Watts distribution of lifetimes for a given parameter set
% using a series aproximation.
%
% INPUT PARAMETERS:
%     tauKWW      = The fitted tauKWW parameter - can be a 1D array (see Indexing)
%     beta        = The fitted beta parameter (0<beta<1) - can be a 1D array (see Indexing)
%     Nterms      = The number of terms to calculate for the series
%     timescale   = The timescale (string) for the rate constants
%     tmin,tmax   = Minimum and maximum times to plot (default: 0,10*tauKWW)
%     plotscale   = Defines the X and Y axes plot scale (e.g. 'linlog', 'loglog', etc.)
%     ShowTauAvg  = Set to 1 to show the value of <tau> on the plot with an arrow
%     Indexing    = Set to 1 to individually index (tauKWW,beta) pairs.
%                   When set to 0 will plot each element of tauKWW versus
%                   each element of beta.
%
% Ricardo J. Fernández-Terán / v1.0 / 26.11.2017

% Initialise defaults
if isempty(tmin)
    tmin = 0.01;
end
if isempty(tmax)
    tmax = tauKWW*10;
end
if isempty(plotscale)
    plotscale = 'linlin';
end
if isempty(ShowTauAvg)
    ShowTauAvg = 0;
end

% Create figure and hold
figure();
hold on

% Calculate and plot, then disable hold
syms k t
L = length(beta);
M = length(tauKWW);
if L==1
    cmap=[0 0 1];
else
    cmap=colormap(jet(L))*0.9;
end
y={};
if indexing==0
    for j=1:M
        for i=1:L
            f=(((-1)^k)/factorial(k))*sin(pi*beta(i)*k)*gamma((beta(i)*k)+1)*(t/tauKWW(j))^((beta(i)*k)+1);
            series = symsum(f,k,0,Nterms);
            rho =(-tauKWW(j)/(pi*t))*series;
            if ShowTauAvg
                tauavg(j) = tauKWW(j)/beta(i)*gamma(1/beta(i));
                name = ['\beta = ' num2str(beta(i)) '; \langle\tau\rangle = ' num2str(tauavg(j),2) ' ' timescale];
                % Add an arrow pointing at <tau>
                t=tauavg(j); rhoatavgtau(i) = double(subs(rho)); clear t; syms t;
                txt = {'\uparrow';'\langle\tau\rangle'};
                text(tauavg(j),rhoatavgtau(i),txt,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',14,'FontWeight','bold','Color',cmap(i,:)*j/M);
            else
                name = ['\beta = ', num2str(beta(i)) '; \tau_{KWW,' j '} = ' num2str(tauKWW(j))];
            end
            fp = fplot(rho,[tmin,tmax],'Color',cmap(i,:)*j/M,'LineWidth',1,'DisplayName',name,'MeshDensity',100);
            y{i,j} = fp.YData;
        end
    end
else
    rhoatavgtau=zeros(L);
    tauavg=zeros(M);
    for i=1:L
        j=i;
        f=(((-1)^k)/factorial(k))*sin(pi*beta(i)*k)*gamma((beta(i)*k)+1)*(t/tauKWW(j))^((beta(i)*k)+1);
        series = symsum(f,k,0,Nterms);
        rho =(-tauKWW(j)/(pi*t))*series;
        switch ShowTauAvg
            case 0
                name = ['\beta = ', num2str(beta(i))];
            case 1
                tauavg(j) = tauKWW(j)/beta(i)*gamma(1/beta(i));
                name = ['\beta = ' num2str(beta(i)) '; \langle\tau\rangle = ' num2str(tauavg(j),2) ' ' timescale];
                % Add an arrow pointing at <tau>
                t=tauavg(j); rhoatavgtau(i) = double(subs(rho)); clear t; syms t;
                txt = {'\uparrow';'\langle\tau\rangle'};
                text(tauavg(j),rhoatavgtau(i),txt,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',14,'FontWeight','bold','Color',cmap(i,:)*j/M);
            case 2
                name = ['\beta = ', num2str(beta(i)) '; \tau_{KWW,' num2str(j) '} = ' num2str(tauKWW(j))];
            case 3
                tauavg(j) = tauKWW(j)/beta(i)*gamma(1/beta(i));
                name = ['\beta = ' num2str(beta(i)) '; \langle\tau\rangle = ' num2str(tauavg(j),2) ' ' timescale];
        end
        fp = fplot(rho,[tmin,tmax],'Color',cmap(i,:)*j/M,'LineWidth',2,'DisplayName',name,'MeshDensity',100);
        y{i,j} = fp.YData;
    end
end
            
% Format plots
hold off
legend('Location','northwest')
legend('boxoff')
legend({},'FontSize',12)
box on
ax=gca;
switch plotscale
    case 'loglog'
        ax.XScale='log';
        ax.YScale='log';
        linlog = 0.001;
    case 'linlin'
        ax.XScale='lin';
        ax.YScale='lin';
        linlog = 0;
    case 'linlog'
        ax.XScale='lin';
        ax.YScale='log';
        linlog = 0.001;
    case 'loglin'
        ax.XScale='log';
        ax.YScale='lin';
        linlog = 0.001;
end
xlim([tmin tmax]);
ylim([linlog max([y{:}])*1.1])
xlabel(ax,['\tau (' timescale ')'],'FontSize',14,'FontWeight','bold');
ylabel(ax,'G_{KWW} (\tau)','FontSize',14,'FontWeight','bold');