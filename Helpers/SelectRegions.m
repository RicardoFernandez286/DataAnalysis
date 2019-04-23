function [RegionPositions,N_regions] = SelectRegions(Axes)
% This function allows to select multiple regions in a graph and output the position of the maximum rectangles.
% The output is a matrix containing the maxima and minima in X and Y of the regions defined by a getline
% commmand.
%
% USAGE:
% Left-click to add points
% Press Return key when done.
% Press backspace to remove last point
%
% OUTPUT:
%     RegionPositions = A matrix with elements of the form [x_min x_max y_min y_max]. Each row defines a region.
%     N_regions       = The number of regions selected by the user.
% 
% Ricardo Fernández-Terán
% v1.0a - 08.06.2018

i=0;
hold on
while 1 % Repeat the loop until Return is pressed
    [xi,yi] = getline(Axes,'closed');
    if length(xi) > 4
        i=i+1;
        x_min   = min(xi);
        x_max   = max(xi);
        y_min   = min(yi);
        y_max   = max(yi);
        RegionPositions(i,:)  = [x_min x_max y_min y_max];
    else
        break
    end
end

N_regions = size(RegionPositions,1);