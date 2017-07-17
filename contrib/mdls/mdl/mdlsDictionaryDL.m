%
% purpose: obtain the description length of a dictionary.
%
% input:
% D ........ Dictionary to be described.
% method ... 1: use ball-n iid model, 2: use laplacian model
%
function [LD,DP]=mdlsDictionaryDL(D,prec,prior,params)
[d,k]=size(D);
% describe length of dictionary
% not needed now
w=sqrt(d);
if ~exist('prior','var')
  prior='Ball';
end
if ~exist('prec','var')
  prec=10; % 10 bits to describe each coefficient
end
if strcmp(prior,'Ball')
  % treat all dictionary coefficients as iid samples from a
  % ball-n distribution, with n=dimension of atoms.
  f=mdlsBallEval(D,d);
  LD=-sum(sum(log2(f))) + d*k*prec;
  DP=D;
elseif strcmp(prior,'Laplacian')
  % describe each atom predictively using the previous neighbor
  % as the prediction value. 
  % to be fair, do it in a zig-zag manner
  DP=mdlsDicPredict(D);
  %
  % ML parameter
  % 
  %L=L(DP|lambda) + L(lambda) + d*k*prec 
  lambda=d*k/sum(sum(abs(DP)));
  %LD=d*k*(-log2(lambda/2)+1)+ d*k/log(2) + 0.5*log(d*k) + d*k*prec;
  fprintf('lambda=%f sqrt(d)=%f\n',lambda,sqrt(d));
  % one option
  %lambda=sqrt(d); % based on ball-n prior on d_ij this is lambda
  LD=d*k*(-log2(lambda/2)+1)+ lambda*sum(sum(abs(DP))) + d*k*prec;
else
  error('Unknown prior specified.');
end
