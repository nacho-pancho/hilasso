%
% function model=mdlsModelSpec(D,alpha,sparsity,lambda,error_metric,coef_metric)
%
% dictionary .. d x k, the codebook or dictionary consisting of k
%               d-dimensional atoms.
% sparsity .... if greater than 0, take only the 'sparsity' absolutly 
%               largest coefficients in the reconstruction of each
%               patch.
% lambda ...... multiplier for the coefficient  penalty term.
% metric ...... select the Lp norm as the metric for the error
%               term, where p=metric.
%
% the struct encapsulates all the information that defines
% a class of sparse models.
%
function model=mdlsModelSpec(width,overlap,dic_size,sparsity,lambda,metric)
  model=struct('width',width,'overlap',overlap,'dic_size',dic_size,'sparsity',sparsity,'lambda',lambda,'metric',metric);
end

