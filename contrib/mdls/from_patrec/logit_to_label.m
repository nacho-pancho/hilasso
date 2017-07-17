%
% INPUT
%
% Y ................ (nxk) logit vector.
%                    each row is of the form [ 0 0 ... 0 1 0 ... 0 ]
%                    where there is a single 1 at the location specificed
%                    by the index of the corresponding label.
%
% OUTPUT: 
%
% y ................ (nx1) each entry corresponds to the index of the
%                    maximum value in each row of Y.
%
function y=logit_to_label(Y)
    y  = max(Y,[],2);
end