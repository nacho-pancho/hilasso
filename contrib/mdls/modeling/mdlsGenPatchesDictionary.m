%
% function D = mdlsGenPatchesDictionary(X,K,sigma,layers,seed)
%
% input:
%
% X ....... patches
% K ....... size of dictionary
% sigma ... std. dev of noise, in 0-255 scale. defaults to 0.
% layers .. how many patches to add to form each initial atom
%
% output:
%
% D ....... the dictionary
%
function D = mdlsGenPatchesDictionary(X,K,sigma,layers,seed)
  if ~exist('sigma','var')
    sigma=0;
  end
  if ~exist('layers','var')
    layers=1;
  end
  if ~exist('seed','var')
    seed=0;
  end
  [n,N]=size(X);
  if seed == 0
      seed = N+n;
  end
  if N == 0
      error('No patches for initializing dictionary');
  end
  rand('twister', seed);
  randn('state', seed);
  D = zeros(n,K);
  bad=1:K;
  while length(bad) > 0      
      for i=1:layers      
          N2 = randperm(N);
          while length(N2) < length(bad)
              N2 = [N2 N2];
          end
          selected = N2(1:length(bad)); 
          D(:,bad) = D(:,bad) + X(:,selected);
      end
      bad = find(sum(D.*D) == 0);
  end
  D = D * (1/layers) + randn(n,K) * sigma;
  for l=1:K
      if norm(D(:,l)) == 0
          error('oops');
      end
    D(:,l) = D(:,l) / norm(D(:,l));
  end

end
