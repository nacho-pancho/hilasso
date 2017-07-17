%
% function model = learn_mcsvm(X,Y)
% 
% Learns a multi-class SVM from samples X and labels Y
%
% X ............ (nxp) n d-dimensional training samples 
% Y ............ (nx1) n training labels, with k different possible values.
% C ............ (1x1) SVM regularization parameter
% debug ........ (bool) whether show debugging or not
%
% OUTPUT
% model ........ learnt model, which is a 1xk cell where each element is 
%                a two-class svm for distinguishing between each class
%                and all the others.
%
function model = learn_mcsvm2(X,y,ktype,kpar,C)
%
% we learn k different SVMs, each one corresponding to a 'one against
% the rest' classifier.
%
    l = sort(unique(y));
    k = length(l);
    model = cell(1,k);
    for i=1:k
        %        fprintf('Learning SVM for class %d ...\n',i);
        fprintf('.');
        li = l(i);
        y2 = 2*(y == li) - 1;
        model{i} = learn_svm(X,y2,ktype,kpar,C);
    end
end
