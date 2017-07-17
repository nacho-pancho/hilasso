function B=mdlsSparsifySlow(A,L) 
  B=zeros(size(A));
  [As,I]=sort(A,'descend');
  for i=1:size(A,2)
      B(I(1:L,i),i)=As(1:L,i);
end