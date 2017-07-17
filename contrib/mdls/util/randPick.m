%%
%% random pick of training samples 
%% useful for bootstrap statistical significance methods.
%%
%% inputs:
%% nsamples ...... total number of samples
%% k ............. number of 
%% npicks ........ if specified, repeat the random pick npicks times
%%
%% outputs:
%% chosen ........ matrix of npicks x k where each i-th row is vector of  k indexes 
%%                 corresponding to the randomly picked samples for the i-th pick.
%% 
function chosen=randPick(nsamples,k,npicks)
if nargin<3
  npicks=1;
end
% 
% a fair pick is obtained as follows:
% nleft=n, kleft=k
% for each i=1:n do
%   toss a binary coin with P(1)=kleft/nleft
%   if 1, 
%      choose i as one of the samples, decrement kleft
%   end
%   decrement nleft
% end
chosen=zeros(npicks,k);
for pick=1:npicks
  nleft=nsamples;
  kleft=k;
  j=1;
  for i=1:nsamples
    if (kleft==0)
      break;
    end
    if rand < (kleft/nleft)
      chosen(pick,j)=i;
      j=j+1;
      kleft=kleft-1;
    end
    nleft=nleft-1;
  end
end
endfunction
