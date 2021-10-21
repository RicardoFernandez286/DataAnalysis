function cov2D = covar2D(Xmat,Ymat)
% Calculate the cross-covariance matrix of matrices Xmat and Ymat
% Xmat and Ymat are observation matrices, where each column is one
% observable (i.e. detector pixel) and each row is one observation (i.e. laser shot)
% Ricardo Fernández-Terán / v1.0a / 2021.08.09

dXmat = (Xmat-mean(Xmat,1))';
dYmat = (Ymat-mean(Ymat,1))';
cov2D = dXmat*dYmat'/(size(Xmat,1)-1);
end