function cov2D = covar2D(Xmat,Ymat)
% Calculate the cross-covariance matrix of matrices Xmat and Ymat
% Xmat and Ymat are observation matrices, where each column is one
% observable and each row is one observation
dXmat = (Xmat-mean(Xmat,1))';
dYmat = (Ymat-mean(Ymat,1))';
cov2D = dXmat*dYmat'/(size(Xmat,1)-1);
end