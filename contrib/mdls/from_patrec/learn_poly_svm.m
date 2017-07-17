%
% function model = learn_poly_svm(X,y,ord,C,debug)
% Learn a polynomial SVM
%
% X ............ (nxp) n d-dimensional training samples 
% y ............ (nx1) n training labels (-1,1)
% ord .......... (1x1) polynomial order
% C ............ (1x1) SVM regularization parameter
% debug ........ (bool) whether show debugging or not
%
function model = learn_poly_svm(X,y,ord,C,debug)
    if ~exist('debug','var')
        debug = false;
    end
    if ~exist('C','var')
        C = 1;
    end
    if ~exist('ord','var')
        ord = 4;
    end
    p = size(X,2);
    ktype = 'poly';
    kpar = ord;
    use2norm = 0;
    model = svm(p,ktype,kpar,C,use2norm,'loqo');
    a0 = [];
    model = svmtrain(model,X,y,a0,debug);
end
