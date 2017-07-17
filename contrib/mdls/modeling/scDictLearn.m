%function D=scDictLearn(X,D0,J,lambda,mode)
%
% Wrapper for using Julien's mexTrainDL as if we were calling the old scDictLearn
%
% 
%   Purpose: learn dictionaries using the L1-regularized model:
%
%   min_{A,D} ||X-DA||_F^2 + lambda ||A||_1
%  
% 
%   Inputs:
%    X ....... Signal to approximate, given as an  n x m matrix where n=dimension of samples, m=number of samples
%    D0 ....... Initial dictionary, given as an n x k matrix  where n=dimension of atoms, k=number of atoms
%    J ....... number of iterations of the algorithm
%    lambda .. Penalty term. Try the magical '0.1' 
%    mode .... Selects minimization problem to be solved:
%     0) min_A ||X-D*A||_2 s.t. ||A_i||_1 <= lambda
%     1) min_A ||A||_1 s.t. ||X_i - D*A_i|| <= lambda
%     2) min_A ||X-D*A|| + lambda*||A||
%
function [D,varargout]=scDictLearn(X,D0,J,lambda,mode)

  if nargin == 0
      help scDictLearn
      D=[];
      return
  end
  params=struct();
  if ~isempty(D0)
      params.D = D0;
  end
  params.iter = J;
  params.lambda = lambda;
  if exist('mode','var')
      params.mode = mode;
  end
  D = mexTrainDL_Memory(X,params);

  if nargout > 1
      max_nz = min(size(X,1)-1,size(D,2));
      varargout{1} = scLasso(X, D, max_nz, lambda, 2);
  end

end
