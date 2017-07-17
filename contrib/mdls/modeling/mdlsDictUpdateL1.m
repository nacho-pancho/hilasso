%
% function D=dicLearnL1step(X,D,A,eta)
%
% purpose:
% one step of  L1 optimization of a Dictionary
% see dicLearnL1 for a complete description.
%
% inputs:
% D ........ previous dictionary, d x k
% X ........ data, d x n
% A .... representation coefs. k x n
% eta ...... quasi-norm auxiliary term. Defaults to 1e-5.
%
% outputs:
% D ........ updated dictionary
%
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function D=mdlsDicUpdateL1(X,D,A,eta)
  if ~exist('eta','var')
    eta=1e-5;
  end
  [m,k]=size(D);
  n=size(X,2);
  E=X-D*A;
  E=1./sqrt(E.*E+eta);
  A=A*(E.*X)';
  for i=1:m
    A2=A.*(ones(k,1)*sqrt(E(i,:)));
    B=A1*A2'; % weighted autocorrelation
    D(i,:)=B\A(:,i);
  end
  % renormalize each column to 1
  for j=1:k
    D(:,j)=D(:,j)/norm(D(:,j));
  end


  
