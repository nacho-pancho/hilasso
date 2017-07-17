%
% function [m M] = mdlsNormalizedMutualInformation(clustering,k,labels)
%
% Input:
%
% clustering ......... cluster labels 
% labels ............. true labels
%
% Output:
% m .................. normalized mutual information
% M .................. confusion matrix (not really).
%
function [m M] = mdlsNormalizedMutualInformation(clustering, labels)
  classes = sort(unique(labels));
  NC = length(classes);
  N = length(clustering);
  c = zeros(1,NC);
  for i=1:NC
      idx = find(clustering == classes(i))';
      c(i) = length(idx);
      for j=1:NC
          M(j,i) = sum( labels(idx) == classes(j) );
      end
  end
  
  C = repmat(c,NC,1);
  
  r = sum(M,2);
  R = repmat(r,1,NC);
  
  id = find(M);
  Mc = M(id).*log(N*M(id)./C(id)./R(id));
  I = sum(Mc(:))/N;
  
  id = find(r);
  H1 = -r(id)'*log(r(id)/N)/N;
  
  id = find(c);
  H2 = -c(id)*log(c(id)'/N)/N;
  
  m = 2*I/(H1+H2);
end
