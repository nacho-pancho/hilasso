%
% function D = NormalizeDictionary(D)
%
<<<<<<< .mine
function D = mdlsNormalizeDic(D,remove_dc)
  if nargin < 2
      remove_dc = false;
  end
=======
function D = mdlsNormalizeDic(D)
  warning('This function is obsolete. use mdlsDictNormalize instead.');
>>>>>>> .r404
  [n,K] = size(D);
  for l = 1:K
      d = D(:,l);
      if remove_dc
          d = d-mean(d);
      end
    tmp = norm(d);    
    if isequal(tmp,0)
        d = randn(n,1);
        d = d-mean(d);
        tmp = norm(d);
    end
    d = d*(1/tmp);
    D(:,l) = d;
  end
end
