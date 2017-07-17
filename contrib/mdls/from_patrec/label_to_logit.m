%
% INPUT
%
% y ................ (nx1) with n arbitrary numeric labels taking k
%                    different values. The labels will be sorted in 
%                    ascending order and an index from 1 to k will be
%                    used to construct the logit output.
%
% OUTPUT: 
%
% Y ................ (nxk) logit-like output: 
%                    each row is of the form [ 0 0 ... 0 1 0 ... 0 ]
%                    where there is a single 1 at the location specificed
%                    by the index of the corresponding label.
%
function Y=label_to_logit(y)
    n  = length(y);
    il = sort(unique(y));    
    k  = length(il);
    Y  = spalloc(n,k,n);
    for i=1:n
        Y(i,find(il==y(i))) = 1;
    end
end