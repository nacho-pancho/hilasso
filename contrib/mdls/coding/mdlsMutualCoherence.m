%
% function varargout = mdlsMutualCoherence(D)
%
% Category: analysis
%
% Description:
%     Computes the various forms of Mutual Coherence for a given dictionary D.
% 
% Input
% D ........ Dictionary to analyze.
% 
% Output
% 
% mc ....... mutual coherence = cosine of minimum angle between any two atoms
% cc ....... cumulative coherence (optional)
%
function varargout = mdlsMutualCoherence(D)

  D=mdlsDictNormalize(D); % coherence is measured on normalized vectors
  K = size(D,2);
  G = abs(D'*D);
  G = G - diag(diag(G));

  varargout{1} = max(abs(G));
  if nargout == 2
    varargout{2} = sum(G,1);
  end

end
