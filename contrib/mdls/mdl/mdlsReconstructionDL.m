%
% purpose: obtain the description length of a dictionary.
%
% input:
% alpha ....... Reconstruction coefficients
% asize ....... Alphabet size (defines range of values)
% prec ........ Quantization step (precision)
% use_ml ...... If zero, use the maximum likelihood estimator 
%               of Lambda instead of the given value.
%               
function LA=mdlsReconstructionDL(alpha,asize,prec,spar,prior,params)
if ~exist('prior','var')
  prior='Laplacian';
end
if ~exist('params','var')
  params=0;
end

[k,n]=size(alpha);
LA=0;
%
% hiperparameter L(sparsity)
%
LA = LA + log2(k); % sparsity takes values in [0,k] -> log2(k) bits
%
% a) nonzero support 
%
% for a given sparsity s the support of the nonzero coefficients
% in each column of alpha is specified by an index with description 
% length of log(k choose s) nats.
% If no sparsity is assumed, there is no index cost.
%
if (spar < k) && (spar>0)
  % Stirling approximation to nchoosek
  indexcost=log(1/sqrt(2*pi))+(k+0.5)*log(k)-(k-spar+0.5)*log(k-spar)-(spar+0.5)*log(spar);
else
  indexcost=0;
end

LA = LA + n*indexcost;

if strcmp(prior,'Laplacian')
  lambda=params(1,1);
  signcost=spar*log(2); % need it in nats for now
  %fprintf('indexcost=%d signcost=%d\n',indexcost,signcost);
  LA = LA + n*signcost;
  %
  % b) PDF evaluated at alpha
  %
  RSA=sum(sum(abs(alpha))); % Residual Sum of Absolute values
  if lambda==0
    lambda=n*spar/RSA;
  end
  fprintf('lambdaML=%f\n',lambda);
  LA = LA - n*spar*log(lambda) + lambda*RSA;
  % parameters: they cannot simply be described with sqrt(n*d)
  % precision!, this is not a 'nice' familiy....
  LA = LA + n*spar*prec; 
elseif strcmp(prior,'GamLap')
  LA=0; % not implemented
else
  error('Unknown prior specified.');  
end