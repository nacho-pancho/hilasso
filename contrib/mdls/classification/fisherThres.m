%%
%% Find the Fisher threshold for the desierd True Positive rate
%%
%% input:
%%
%% fld .... Projection vector
%% X ...... class 1 data, each column corresponds to the features of a single sample
%% Y ...... class 2 data, idem.
%% tp ..... desired true positive rate, 0 < tp <= 1
%%
%% output:
%%
%% tau .... threshold to be used
%%
function [tau,fp] = fisherThres(fld,X,Y,tp)

[Dim Nx]=size(X);
[Dim Ny]=size(Y);

% value depends only on Y

projY=fld'*Y;
projX=fld'*X;
[sortedY,i]=sort(projY);
desired_n=Ny*(1-tp);
if floor(desired_n)==desired_n % is an integer
  tau = sortedY(desired_n);
else
  above = ceil(desired_n);
  below = above-1;
  r   = desired_n-below;
  tau = sortedY(below)+r*(sortedY(above)-sortedY(below)); % linear interp.
end
fp = length(find(projX > tau))/Nx;
endfunction