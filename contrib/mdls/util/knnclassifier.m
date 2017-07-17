function [label,c,kNNlabel] = knnclassifier(X,D,k)
%
% Function that classifies acording to de KNN
%
%   Input:
%
%   D :  is a cell containing the training samples pre class
%   X :  n x m matrix, of the data to be analyzed (m n-dim vectors)
%   k :  is a vector that contians the number of neighbors to consider
%
%
%   Output:
%
%   c    : 

% number of classes
C = length(D);
% number of values of k
K = length(k);

Dr = [];
len = zeros(1,C);
for i=1:C
    Dr = [Dr D{i}];
    len(i) = length(D{i}); 
end

cumlen = [0 cumsum(len)];

kmax = max(k);

% Compute kNN
idx = knnsearch(X',Dr',kmax)';

c = zeros(K,C);

for j=1:K
    
    % Do the voting
    kNNlabel = zeros(size(idx(1:k(j),:)));
    for i=1:C
        kNNlabel( idx(1:k(j),:) > cumlen(i) & idx(1:k(j),:) <=cumlen(i+1)) = i;
    end

    % assign class
    label = mode(kNNlabel);
    
    % cumpute number of elements per class
    c(j,:) = hist(label,1:C);

end

    
    
    