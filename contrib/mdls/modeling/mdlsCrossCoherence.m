%
% Measures cross-coherence betweeen dictionaries. Given C dictionaries,
% the function returns two matrices of size (C+1)x(C+1).
% In the first one, the element c_ij  contains the average
% cross-coherence between  dictionaries i and j.
% In the second, the element c_ij contains the maximum atom-wise cross
% coherence between dictionaries i and j.
% The diagonal elements correspond to the mutual coherence of each
% dictionary.
%
% The elements (i,end) and (j,end) of both matrices contain the row and
% column wise averages respectively, without taking into account the diagonals.
% Likewise, the element (end,end) contains the average of the
% off-diagonal element of the whole matrix.
%
% Input:
% D ......... cell where each element is a dictionary
%
% Output:
% mxc ....... mean cross-coherence matrix
% Mxc ....... maximum cross-coherence matrix
%
function [mxc,Mxc] = mdlsCrossCoherence(D)
  C=length(D);
  mxc = zeros(C+1,C+1);
  for i=1:C
    for j=1:C
      M=D{i}'*D{j};
      mxc(i,j) = mean(abs(M(:)));
      Mxc(i,j) = max(abs(M(:)));
    end
  end
  mxc(:,end) = (sum(mxc,2)-diag(mxc))/(C-1);
  Mxc(:,end) = (sum(Mxc,2)-diag(Mxc))/(C-1);
  mxc(end,:) = mxc(:,end)';
  Mxc(end,:) = Mxc(:,end)';
  mxc(end,end)=mean(mxc(1:(end-1),end));
  Mxc(end,end)=mean(Mxc(1:(end-1),end));
end
