%
% function model = learn_linear_svm(X,y,C,debug)
% Learn a Linear SVM
%
% X ............ (nxp) n d-dimensional training samples 
% y ............ (nx1) n training labels (-1,1)
% ktype ........ (char) kernel type. See kernel from svm package.
% kpar ..........(1x?) kernel parameter(s) corresponding to the chosen
%                      kernel type.
% C ............ (1x1) SVM regularization parameter
% debug ........ (bool) whether show debugging or not
%
% OUTPUT
% model ........ learnt model.
%
function model = learn_svm(X,y,ktype,kpar,C,debug)
    if ~exist('debug','var')
        debug = false;
    end
    if ~exist('ktype','var')
        ktype = 'linear';
    end
    if ~exist('kpar','var')
        kpar = [];
    end
    if ~exist('C','var')
        C = 1;
    end
    p = size(X,2);
    use2norm = 0;
    model = svm(p,ktype,kpar,C,use2norm,'loqo');
    a0 = [];
    model.qpsize = 100;
    model = svmtrain(model,X,y,a0,debug);
end
