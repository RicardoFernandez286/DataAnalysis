function plotAngleHist(simData,WhichTraject,HistPlot)
%% READ from simData
version				= simData.version;
trajectInfo			= simData.trajectInfo;
trajectData_2D 		= simData.trajectData_2D;
trajectData_3D 		= simData.trajectData_3D;
Nmolecules			= simData.Nmolecules;
Nsamples			= simData.Nsamples;
Radius_LJ_Re		= simData.Radius_LJ_Re;
Radius_LJ_CN		= simData.Radius_LJ_CN;
Rbox				= simData.Rbox;

while WhichTraject == 0
    WhichTraject = round(Nsamples*rand);
end
            
% Isotope shifts
iso13_shift = 48;
iso18_shift = 92;

% Dipole moment magnitude (!)
mu          = mean(sqrt(sum((trajectData_2D(trajectData_2D(:,6) == -iso13_shift,3:5)).^2,2)));


switch HistPlot
        case 'Single'
            trajData    = squeeze(trajectData_3D(:,:,WhichTraject));
            trajData    = trajData(:,3:5);
        case 'All'
            trajData    = trajectData_2D(:,3:5);
end


angles  = acos(trajData(:,3)./sqrt(sum(trajData.^2,2)))*180/pi;

fh3     = figure(3);
clf(fh3)
histogram(angles,100)
title('Angular distribution histogram')
xlabel('Inclination angle, \theta (\circ)','FontWeight','bold','FontSize',14)
ylabel('Counts','FontWeight','bold','FontSize',14)
xlim([0 90]);

end