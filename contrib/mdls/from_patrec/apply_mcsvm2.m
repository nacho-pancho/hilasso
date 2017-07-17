%
% function y = apply_mcsvm2(model,X)
% 
% Apples a multi-class SVM to samples X
%
% X ............ (nxp) n d-dimensional training samples 
% model ........ model, which is a 1xk cell where each element is 
%                a two-class svm for distinguishing between each class
%                and all the others.
% debug ........ (bool) whether show debugging or not
%
% OUTPUT
%
% y ............ (nx1) output labels.
%
% PENDING: could be improved with some confidence level by
% comparing the different classifications instead of keeping the class
% that is farthest form the margin
%
function y = apply_mcsvm2(model, X)
 k = length(model);
 n = size(X,1);
 Y = zeros(n,k);
 for i=1:k
     [dum,Y(:,i)] = svmfwd(model{i},X);
 end
 [dum,y] = max(Y,[],2);
end