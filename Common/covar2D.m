function cov2D = covar2D(Xmat,Ymat)
% Calculate the cross-covariance matrix of matrices Xmat and Ymat
% Xmat and Ymat are observation matrices, where each column is one
% observable and each row is one observation

Rmat        = corr(Xmat,Ymat);
Sigma_mat   = std(Xmat,0,1)'*std(Ymat,0,1);
cov2D       = Rmat.*(Sigma_mat);

end