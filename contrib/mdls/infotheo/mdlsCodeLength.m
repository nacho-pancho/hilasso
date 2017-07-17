%
% Given a set of observed frequencies f and given probabilities p
% for a set of symbols, this computes the total codelength
% obtained with Shannon coding.
%
% Warning: if p contains zeros where f has not, then
% this will give nice (and expected) inf codelength
% unless clip = true, in which case the 'unexpected' symbols
% will be encoded losslessly as the closest symbol with positive probability.
% Assuming unimodal probabilities.
%
function L=mdlsCodeLength(f,p,clip)
    if nargin < 3
        clip=false;
    end
    if sum(f(p==0)) > 0
        warning(['Nonzero frequencies for symbols with zero ' ...
                 'probabilities!.']);
        if (clip) 
            ppi = find(p>0);
            maxpos = max(ppi);
            minpos = min(ppi);
            pfi = find(f>0);
            right_tail = pfi(pfi>maxpos);
            left_tail = pfi(pfi<minpos);
            f(maxpos)=f(maxpos)+sum(f(right_tail));
            f(right_tail)=0;
            f(minpos)=f(minpos)+sum(f(left_tail));
            f(left_tail)=0;
        else            
            L=inf;
            return;
        end
    end
    nzi=find(p~=0);
    L=-sum(f(nzi).*log2(p(nzi)));
end