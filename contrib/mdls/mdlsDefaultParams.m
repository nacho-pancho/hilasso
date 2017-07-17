%
% function params = DefaultParams()
%
% Name : DefaultParams
%
% Category: utility
%
% Output:
% params ....... (struct) model parameters with the following members:
% MODEL PARAMETERS:
%                patchWidth ........ Width of patches
%                patchOverlap ........... Patch overlapping
%                varError .......... Variance of error term for Gaussian priors.
%                betaError ......... MOL beta parameter for error term.
%                kappaError ........ MOL kappa parameter for error term.
%	         betaCoef .......... MOL parameter for reconstruction
%                                    coefficient regularization term.
%	         kappaCoef ......... MOL kappa parameter for  reconstruction.
%                lambdaCoef ........ Laplacian parameter for rec. coef.
%                thetaCoef ......... Bernoulli param for rec. coeff.
%                mu ................ Penalty term for dictionary 
%                                    regularization term.
%                eta ............... Penalty term for atom norms. Use '0' for
%                                    original cost function (and enable
%                                    forceNormalization), or use a large value
%                                    (say greater than 50).
%                lambda ............ L1 penalty. If 0 (default) it is computed
%                                    from eta, mu and varError
%                modelingIters ..... Total iterations for model training
%                modelingTolerance . Threshold to stop iterations in modeling
%
% CODING:
%                codingIters ....... Iterations for coding step
%                codingMode ........ Sparse coding mode: 
%                                    0: constrained error, 
%                                    1: constrained regularization term, 
%                                    2:lagrangian (DEFAULT)
%                codingMethod ...... Method for sparse coding. Can be
%                                    - Lasso : Lasso solution using LARS
%                                    - MOL   : Mixture of Laplacians
%                codingSubMethod ... Only MOD: method used for the weighted L1
%                                    minimization performed at each iteration.
%                                    Default 'Lasso'. Can also be 'IS'.
%                codingTolerance ... Threshold to stop coding step.
%                codingSparsity .... Desired sparsity in coding output
%
% UPDATE:
%                updateMethod ...... Method for dictionary update, can be:
%                                    - MOCOD: fixed point MOCOD-like update 
%                                    - ND : Newton MEthod
%                                    - SD : Plain Steepest Descent
%                                    - CG : Conjugate Gradient
							     %                updateStepsizeRule  Stepsize rule: ARO
%                updateIters ....... Iterations for dic update step
%                updateTolerance ... Threshold to stop dic. update step.
%                reviveAtoms ....... Reinitialize atoms that are unused.
% MISCELANEOUS
%                debug ............. Show debug info and images
%                forceNormalization  If different from 0, force a 
%                                    normalization step every 
%                                    'forceNormalization' update steps.
%
function params = mdlsDefaultParams()
params = struct(...
  'patchWidth',8,...
  'patchOverlap',4,...
  'costFunction','orig',...
  'varError',.0001,...  % for normalized input data this is 1% error
  'betaError',10,... 
  'kappaError',2,...
  'betaCoef',.1,...
  'kappaCoef',3,...
  'lambdaCoef',(3-1)/.1, ... % mean mol = (kappa-1)/beta
  'thetaCoef',0.1,... % sparsity = 10%
  'mu',0.1,...
  'eta',50,...
  'lambda',0,...
  'modelingIters',10,...
  'modelingTolerance',1e-3,...
  'codingMethod','Lasso',...
  'codingMode',2,...
  'codingSubMethod','Lasso',...
  'codingIters',10,...
  'codingTolerance',1e-3,...
  'codingSparsity',5,...
  'codingIncremental',false,...
  'updateMethod','MOL',...
  'updateIters',10,...
  'updateTolerance',1e-2,...
  'updateStep',1e-2,...
  'updateStepsizeRule','ARO', ...
  'updateIncremental',true,...
  'reviveAtoms',false,...
  'debug',1,...
  'outdir','../output',...
  'forceNormalization',0);
end

