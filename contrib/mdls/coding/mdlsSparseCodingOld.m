
% function A = mdlsSparseCoding(X,D,A0,params)
%
% Name: mdlsSparseCoding
%
% Category: core function
%
% Description: wrapper for LASSO amd MOL under a single interface.
%
% Input:
% X ........ data (n x N)
% D ........ dictionary (n x K)
% A0 ....... initial coefficients (K x N). If a scalar of value a is passed
%            instead, A0 will be computed using the penalized least squares
%            solution with penalty a.
% params ... see PDSparseModeling for details
%
% Output:
% A ........ updated reconstruction coefficients
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%
function A = mdlsSparseCoding(X,D,A,params)
  %
  % the variance of the error is the parameter of the fitting term which
  % to simplify the problem gets factorized out, which means it multiplies
  % the value of Lambda
  %
  sigma2 = params.varError;
  b = params.betaCoef;
  k = params.kappaCoef;
  lambda = params.lambda;
  theta = params.thetaCoef;

  %fprintf('Equivalent L1 LAMBDA=%f\n',lambda);
  method = params.codingMethod;

  %
  % if one wants 'equivalent' parameterizations between MOL and Laplacian
  % the relationship between both is the mean value that the Gamma
  % prior 
  % 
  % Gamma(theta;kappa,beta)
  % 
  % gives for the Exponential parameter theta
  %
  % (Exp(x;theta)=theta*exp(-theta*x))
  %
  % which is theta = kappa/beta
  %
  % so the equivalent L1 penalty is : lambda = 2*sigma^2*theta = 2*sigma^2*kappa/beta
  % 
  % on the other hand, for fixed lambda and sigma^2 we have theta =
  % lambda/(2*sigma^2). Using theta = kappa/beta => beta = kappa/theta
  % we have the equivalent MOL 
  %
  % beta = kappa/(lambda/2*sigma^2) = 2*sigma^2*kappa/lambda
  %
  % This is implicitly set by setting lambda = 0 or beta = 0
  %
  if (lambda == 0) && (b == 0) && (theta == 0)
      error(['Must provide valid values for at least one of, params.lambdaCoef, params.thetaCoef ' ...
             'or params.betaCoef']);
  end

  if b == 0
      if lambda ~= 0
          b = 2*sigma2*k/lambda;
      else
          b = k/theta;
      end
  end  
  sparsity = params.codingSparsity;
  N = size(X,2);
  K = size(D,2);
  Ns = ceil(N/10);
  nd = 0;
  incremental = params.codingIncremental;

  if isequal(method,'IS')
    %
    % -------- ITERATED SHRINKAGE
    %
    c = 1.1*norm(D'*D);
    if lambda == 0
        if theta ~= 0 
            lambda = 2 * sigma2 * theta;
        else
            lambda = 2 * sigma2 * k / b; 
        end
    end
    %c = 1;
    %fprintf('IS constant c=%f\n',c);
    A = LeastSquaresSparseCoding(X,D,1e-1,2);
    A = scIteratedShrinkage(X,D,A,c,lambda,params.codingIters,params.codingTolerance);
  elseif isequal(method,'Lasso')
    %
    % -------- LASSO
    %
    if lambda == 0
        if theta ~= 0 
            lambda = 2 * sigma2 * theta;
        else
            lambda = 2 * sigma2 * k / b; 
        end
    end
    A = scLasso(X, D, sparsity, lambda, params.codingMode);
  else % end LASSO
    %
    % -------- MOL
    %
    % MOL's lambda is lambda=sigma^2*(k+1)
    if (params.codingMode == 2)
        lambda = 2 * sigma2 * (k+1);
    end
    A = scMOL(X, D, sparsity, lambda, b, params.codingMode, params.codingIters, params.codingTolerance);
  end % end if MOL
  for i = 1:nd
    fprintf('\b');
  end
end % function
