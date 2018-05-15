function [X,Y,Z] = xyz2mat(xyz)
% USAGE: [X,Y,Z] = xyz2mat(xyz)
% INPUTS:
%     xyz = (MxN)x3 double
% OUTPUTS:
%     X   = Mx1 list of the unique elements in the first column of XYZ
%     Y   = Nx1 list of the unique elements in the second column of XYZ
%     Z   = An MxN matrix, sorted according to X and Y.
%     
% Ricardo Fernández-Terán, v1.0 / 23.03.2018

X       = unique(xyz(:,1));
Y       = unique(xyz(:,2));
Z       = zeros(length(X),length(Y));
z_list  = xyz(:,3);

for i=1:length(X)
    for j=1:length(Y)
        element = (xyz(:,1)==X(i))&(xyz(:,2)==Y(j));
        if cumsum(element) == 0
            Z(i,j) = NaN;
        else
            Z(i,j)  = z_list(element);
        end
    end
end