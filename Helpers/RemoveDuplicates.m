function [Xnew,Ynew] = RemoveDuplicates(X,Y)
% This function will average the duplicate rows from a matrix Y, considering the unique rows of an index matrix X
% Ricardo Fernández-Terán, v1.0 / 2019.10.22

[Xnew,~,NU_id]  = unique(X);
Nunique         = length(Xnew);
Ynew            = zeros(Nunique,size(Y,2));

for i=1:length(Xnew)
    Ynew(i,:)   = mean(Y(NU_id == i,:),1);
end
