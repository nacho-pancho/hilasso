% 
% function [D,alpha,DC]=mdlsCreateModel(X,K,iters,lambda,metric,annealing,outdir)
%
% purpose: create a sparse model for a given set of patches.
%
% input:
% X ........... dxn matrix representing n d-dimensional input vectors, usually unrolled image patches.
% K ........... if scalar, specifices number of atoms in dictionary, else, it is taken as the initial dictioanry itself.
% iters ....... algorithm iterations.
% lambda ...... penalty term used in training
% metric ...... p=metric uses the p-norm to measure reconstruction error.
% annealing ... add sigma=annealing amount of Gaussian noise to initial patches in dictioanry.
% outdir ...... is specified, diagnostic information will be written to this directory.
%
% output:
% D ........... learned dictioanry
% alpha ....... obtained corresponding reconstruction coefficients
% DC .......... if requested as third variable, the DC will be removed from each i-th patch in X and
%               returned as the i-th element in this 1xn vector.
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function [D,alpha,DC]=mdlsCreateModel(X,K,iters,lambda,metric,annealing,outdir)
mdlsConstants; % load various project constants
if ~exist('annealing','var')
  annealing=0;
end
[d,n]=size(X);
if nargout > 2
  % remove DC
  DC=sum(X)/d;
  X=X-ones(d,1)*DC;
end
sigma=annealing;

% the second parameter can be the size of the dictionary, or a given initial dictionary
if size(K,1)*size(K,2) > 1
%fprintf('Initializing dictionary with perturbation sigma=%f.\n',sigma);
  Dini=mdlsDicIni(X,K,sigma,0);
else
  Dini=K;
  K=size(Dini,2);
end

%fprintf('Training dictionary using L%d metric\n',metric);
if metric==2
  fprintf('Creating L2 model using mexDL ... ');
  [D,alpha]=mexDL(X,Dini,iters,lambda,2); % last param is mode=2
  fprintf('\n');
else
  fprintf('Creating L1 model using mdlsDicLearnL1 ... ');
  method=11;
  [D,alpha]=mdlsDicLearnL1(X,Dini,iters,lambda,1e-2,method.outdir);
  fprintf('\n');
end
