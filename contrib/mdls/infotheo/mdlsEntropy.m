%
% Returns the entropy of a random variable given the probabilities
% of its values.
%
function H=mdlsEntropy(p)
    p=p(p~=0); % 0log0 = 0
    H=-sum(p.*log2(p));
end