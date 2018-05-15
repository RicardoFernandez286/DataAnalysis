function LoadAndPlotStreak(data)
% Get time and wavelength vectors
time = data(:,1);
wavelength = transpose(data(1,:));
% Remove first row and first column to get only Z data
data(1,:)=[]; data(:,1)=[];
time(1)=[];wavelength(1)=[];
% Make contour plot
contourf(time, wavelength, data, 10,'LineStyle','-','LineColor','flat');
% Invert the view (optional)
view(-90,90)
set(gca, 'ydir', 'reverse');