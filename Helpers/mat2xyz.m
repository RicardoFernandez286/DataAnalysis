function [X,Y,Z] = mat2xyz(M)
% USAGE: [X,Y,Z] = mat2xyz(mat)
% INPUTS:
%     mat = MxN matrix
% OUTPUTS:
%     X   = List of X indices
%     Y   = List of Y indices
%     Z   = Matrix, sorted out according to the indices in X and Y
%     
% Ricardo Fernández-Terán, v1.0 / 17.07.2019

[X,Y] = ind2sub(size(M), 1:numel(M));
Z = M(:);