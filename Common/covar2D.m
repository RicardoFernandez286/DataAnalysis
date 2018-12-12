function cov2D = covar2D(Xmat,Ymat)
% Calculate the cross-covariance matrix of matrices Xmat and Ymat
% Xmat and Ymat are observation matrices, where each column is one
% observable and each row is one observation

Rmat        = corr(Xmat,Ymat);

[Nrow Ncol] = size(Xmat);

Sigma_mat   = zeros(Ncol,Ncol);
for i=1:Ncol
    Sigma_mat(i,:) = std(Xmat,2)
end


cov
end