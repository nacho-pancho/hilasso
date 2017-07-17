%
% function D = mdlsDictNormalize(D)
%
function D = mdlsDictNormalize(D)
  [n,K] = size(D);
  for l = 1:K
    tmp = norm(D(:,l));
    if isequal(tmp,0)
      D(:,l) = randn(n,1);
      tmp = norm(D(:,l));
    end
    D(:,l) = D(:,l)/tmp;
  end
end
