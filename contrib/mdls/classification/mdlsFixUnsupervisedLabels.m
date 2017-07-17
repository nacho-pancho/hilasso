%
% when performing unsupervised classification, exchange labels so 
% that they correspond to the (presumably) right classes.
% That is, if we get a confusion matrix such as
%
% 1 0 0
% 0 0 1
% 0 1 0
%
% it is likely that the labels are wrong, and that lables of clases
% 2 and 3 should be swapped.
%
function fixed_labels = mdlsFixUnsupervisedLabels(classes,conf_matrix,labels)
  %
  % first build a permutation matrix out of the conf. matrix
  %
  if size(labels,2) == 1
      labels=labels';
  end
  N = size(labels,2);
  conf_matrix = conf_matrix(1:(end-1),1:(end-1));
  C = size(conf_matrix,1);
  [maxhit,mi]  = max(conf_matrix);
  fixed_labels = zeros(size(labels));
  for i=1:N
      ci = find( classes==labels(i) );
      fi = mi( ci );
      fixed_labels(i) = classes(fi);
  end
end