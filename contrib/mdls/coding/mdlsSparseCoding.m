%
% wrapper for various types of sparse coding algorithms.
% The type is selected by the value in params.reg_type.
% The following regularization methods are implemented:
%
% l0 .... l_0 norm regularization. 
%
%         min_A ||X-D*A||_2 s.t. ||A||_0 <= lambda
%         
%         Only mode 0 is implemented
%         here, using Orthogonal Matching Pursuit.
%         
%
% l1 .... l_1 norm regularization. All modes are available.
%
% l2 .... l_2 norm regularization. Only mode 2 for now,
%         but it is straightforward to implement modes 0 and 1.
%
% wl2 .... weighted l_2 norm regularization. Only mode 2 for now,
%         but it is straightforward to implement modes 0 and 1.
%
% wl1 ... weighted l_1 regularization. Here reg(A) = ||W*A||_1
%         where W is a diagonal, positive definite matrix.
%         PARAMETERS: params.theta_row is a kx1 vector with weights for
%         each row.
% 
% rwl2 . this is actually a fast approximation to the l1 solution,
%         specially useful when the solution sought is dense (otherwise
%        it may actually be significantly slower than l1).
%        Only mode 2 is implemented for now. The algorithm approximates
%
%         min_A 0.5*||X-D*A||_2^2 + \lambda ||A||_1 
%
%       using the following iteration:
%
%         A(t+1) = min_A ||X-D*W(t)A||_2^2 + \lambda ||W(t)A||_1 
%
%       where W(t)[k] = 1/(|A(t)[k]| + eps)
%
% moe .. this is the regularizar that results from the Mixtures of
%        Exponentials universal model for sparse coefficients. The
%        regularizer has the following form:
%
%        reg(A) = (kappa+1)\sum{ \log( |A[k]| + \beta) }
%
%        this is a non-convex regularizer whose solution is approximated
%        using an iterative local linear approximation (LLA). 
%        
%        see help scMOL for a detailed explanation.
%
%        PARAMETERS: params.kappa, params.beta
%
% joe .. this derives from the Jeffreys mixture of Exponentials universal
%        model. In this case we have
%     
%        reg(A) = ...
%
%        see help scJOE for a detailed explanation.
%
%        PARAMETERS: params.lambda_min, params.lambda_max
%
% INPUT
%
% X ....... matrix of data samples ordered as coloumns.
% D ....... dictionary of atoms.
% params .. model parameters. 
%   -reg_type ..... choice of the regularizer function.
%   -reg_mode ..... regularization mode: %
%     0) min_A ||X-D*A||_2 s.t. reg(A_i) <= lambda
%     1) min_A reg(A) s.t. ||X_i - D*A_i|| <= lambda
%     2) min_A ||x-D*A|| + lambda*reg(A)
%   -positive ..... if implemented for the selected type, restricts
%                  all coefficients in A to be positive.
%   -some regularizers require extra parameters which
%    have to be provided by the user. See description above for details.
%
% verbose . numeric verbosity level. 0 means nothing.
%
%
%
%
% OUTPUT
%
% A ........ sparse coding coefficients.
%
function A = mdlsSparseCoding(X,D,params,verbose)
  if ~exist('params','var')
      params = mdlsDefaultModelParams();
  end
  if ~exist('verbose','var')
      verbose = 0;
  end
  maxL = min(size(D,1),size(D,2));
  if (params.L <= 0) || (params.L > maxL)
      params.L = maxL;
  end
  switch (params.reg_type)
    case 'l0'
      A = scOMP(X,D,...
                params.L, ...
                params.lambda);
    case 'l1'
      A = scLasso(X,D,...
                  params.L, ...
                  params.lambda, ...
                  params.reg_mode, ...
                  params.positive);
    case 'wl1'      
      A = scRowLasso(X,D,...
                  params.L, ...
                  params.lambda, ...
                  params.theta_row,...
                  params.reg_mode);
    case 'moe'
      A = scMOL(X,D,params.L, ...
                params.lambda, ...
                params.kappa,...
                params.beta,...
                params.reg_mode,...
                params.lla_iter);
    case 'joe'
      A = scJOE(X,D,params.L,...
                params.lambda,...
                params.lambda_min, ...
                params.lambda_max,...
                params.reg_mode,...
                params.lla_iter);

    case 'wl2' 
      %
      % this applies a weighted L2, where the weights
      % are estimated as the empirical variance of the previous
      % run
      DTX = D'*X;
      DTD = D'*D;
      A = inv(DTD+0.5*params.lambda*eye(size(D,2)))*DTX;
      for i=1:params.iter
          w = 0.5*params.lambda./(mean(A.^2,2)+0.1);
          A = inv(DTD+diag(w))*DTX;
      end
    case 'l2' 
      %
      % this is plain old ridge regression
      %
      A = inv(D'*D+0.5*params.lambda*eye(size(D,2)))*D'*X;
    case 'rwl2'
      %
      % this approximates the L1 solution as a series of
      % weighted L2 solutions. This one is better to implement in C
      % due to the loops, but I put it here for now.
      DTX = D'*X;
      DTD = D'*D;
      A = inv(DTD+0.5*params.lambda*eye(size(D,2)))*DTX;
      % unfortunatelly the for kills what would be quite fast...
      for i=1:size(A,2)
          A(:,i) = inv(DTD + params.lambda*diag(0.5./(abs(A(:,i))+eps)))*DTX;
      end
  end
end
