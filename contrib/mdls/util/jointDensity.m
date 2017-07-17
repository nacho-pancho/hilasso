%%
%% compute an empirical estimation of the joint density
%% between random samples x_i and y_i in x and y respectively
%%
function [xy,fxy]=jointDensity(x,y,stepx,stepy)

  %
  % determine ranges
  %
  my=min(y); My=max(y);
  mx=min(x); Mx=max(x);
  %
  % determine sampling step
  %
  Nx=length(x); Ny = length(y);
  if (Nx != Ny)
    error("vectors must have the same length.We have %d <> %d\n",Nx,Ny);
  end
  if nargin < 3
    stepx=1;
    stepy=1;
  end

endfunction