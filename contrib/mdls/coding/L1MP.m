function alpha=L1MP(X,D,s)
[d,k]=size(D);
alpha=zeros(k,1);
res=X;
aux=ones(1,k);
for is=1:s
  %
  % obtain a row of values with the L1 norm of the difference
  % between X and each atom.
  %
  L1d=sum(abs(D-X*aux));
  %
  % select the one that's closer in L1 norm
  %
  [m,im]=min(L1d);
  %
  % add it to solution ... but how??? projection??
  %
  
end
