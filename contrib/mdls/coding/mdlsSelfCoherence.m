%
% function A=selfCoherence(D)
%
% purpose:
% study the structure of a dictionary through the cosine (angle)
% distance between the atoms. The minimum value of the matrix
% returned corresponds to the mutual coherence of the dictionary
% as defined in the sparstiy literature.
% Note: This is simply D'*D.
%
% input:
% D .......... Dictionary to be studied
%
% output:
% A .......... k x k matrix where each column j contains the cosine of the 
% angle of the j-th atom with the other k atoms (including
% himself). 
%
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function A=mdlsSelfCoherence(D)
  A1=abs(D'*D);
  % remove the diagonal zeros
  k=size(A1,1);
  A=zeros(k-1,k);
  for kk=1:k
    if kk > 1
      A(1:(kk-1),kk)=A1(1:(kk-1),  kk);
    end
    if kk < k
      A(kk:end,  kk)=A1((kk+1):end,kk);
    end
  end
  %ki=ones(k,1)*(1:k);
  %A=A( (ki-ki')~=0 );
end
