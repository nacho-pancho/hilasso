%
% function model = learn_linear_svm(X,y,C,debug)
% Learn a Linear SVM
%
% X ............ (nxp) n d-dimensional training samples 
% y ............ (nx1) n training labels (-1,1)
% C ............ (1x1) SVM regularization parameter
% debug ........ (bool) whether show debugging or not
%
% OUTPUT
% model ........ learnt model.
%
function model = learn_linear_svm(X,y,C,debug)
    if ~exist('debug','var')
        debug = false;
    end
    if ~exist('C','var')
        C = 1;
    end
    p = size(X,2);
    ktype = 'linear';
    kpar = [];
    use2norm = 0;
    model = svm(p,ktype,kpar,C,use2norm,'loqo');
    a0 = [];
    model.qpsize = 100;
    model = svmtrain(model,X,y,a0,debug);
end
