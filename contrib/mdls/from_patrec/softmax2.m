%
% safeguarded version to avoid NaNs
%
function S=softmax2(A)
    A = A - repmat(mean(A),size(A,1),1);
    S = exp(A) ./ repmat( sum(exp(A)+eps), size(A,1), 1);
end