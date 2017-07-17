%
% function D = NormalizeDictionary(D)
%
function D = mdlsNormalize(X)
  [M,N] = size(X);
  aux = sqrt(sum(X.*X));
  aux2 = find(aux==0);
  if ~isempty(aux2)
      aux(aux2) = 1;
  end
  X = X.*(repmat(1./aux,M,1));
end
