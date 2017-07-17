%
% purpose: obtain the description length of the DC terms
% of a set of patches.
%
% input:
% dc ........ DC terms to be described
% asize ..... alphabet size
% prec ...... decimal precision in bits, usually log2(dim)
% prior ..... prior to use. 'Uniform', ...
% params .... parameters for the prior.
%
function LDC=mdlsDCDL(dc,prec,prior,params)
[m,n]=size(dc);
if ~exist('prior','var')
  prior='Uniform';
end
if ~exist('prec','var')
  prec=8; % 8 bits to describe each coefficient
end
if ~exist('params','var')
  params=0;
end

if strcmp(prior,'Uniform')
  % continuous uniform distribution in [0,255] discretized in
  % prec steps.
  asize=params(1,1);
  LDC=n*(log(asize)-log(prec));  
else
  error('Prior not implemented.');
end
