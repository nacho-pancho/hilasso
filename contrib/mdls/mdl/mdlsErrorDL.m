%
% purpose: obtain the description length of a dictionary.
%
% input:
% E ........... Error coefficients
% asize ....... Alphabet size (defines range of values)
% prior ....... Prior: 'Gaussian', 'Laplacian','GL' are valid.
%               Defaults to 'Gaussian'.
% params ...... Prior parameters.
%               If zero or empty the maximum likelihood
%               estimations will be used (and encoded).
%               
function LE=mdlsErrorDL(E,precision,prior,params)
[d,n]=size(E);
LE = 0;
if ~exist('prior','var')
  prior='Gaussian';
end
if ~exist('params','var')
  params=0;
end
if strcmp(prior,'Gaussian')
  if isempty(params)
    sigma2=0;
  else
    sigma2=params(1,1);
  end
  %
  % Gaussian prior on error term:
  %
  % n*(d/2)*log 2*pi*sigma2 - 1/2 * 1/sigma^2 * ||E||^2_2
  %
  % taking sigma2 to be the ML estimator
  %
  % sigma2 = 1/(nm) * ||E||^2_2
  %
  % we describe the ML estimator with a precision of 1/sqrt(n),
  % and assuming the maximum possible range of [0,(255/2)^2]
  % yielding a description length of 1/2 log n + 2 log asize/2
  %
  % we obtain
  %
  RSS=sum(sum(E.*E));
  %
  % sigma2=0 indicates the use of a ML estimator
  % (which also has to be described.)
  if sigma2==0
    sigma2=RSS/(d*n);
    % description cost of sigma2.
    % additional cost for describing parameter quantized to
    % sqrt(d*n) precision:
    LE= LE + 0.5*log(d*n); 
    fprintf('sigma_E=%f\n',sqrt(sigma2));
  end
  % PDF
  LE = LE + 0.5*d*n*log (2*pi*sigma2) + 0.5*RSS/sigma2;
elseif strcmp(prior,'Laplacian')
  if params == 0
    lambda=mdlsLaplacianFit(unroll(E));
    fprintf('lambda_E=%f\n',lambda);
    % additional cost for describing parameter quantized to
    % sqrt(d*n) precision:
    LE = LE + 0.5*log(d*n);
  end
  LE = LE - d*n*log(lambda/2) + lambda*sum(sum(abs(E)));
elseif strcmp(prior,'GamLap')
  LE=0; % not implemented
else
  error('Unknown prior specified.');
  LE=-1;
end
%
% this is due to the precision in quantizing the continuous PDF
%
LE = LE + d*n*precision;
