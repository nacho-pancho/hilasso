%
% function DP=mdlsDicPredict(D)
%
% purpose:
% 
% Represent a dictioanry predictively using a simple sequential first order
% linear predictor.
%
% inputs:
%
% D ......... d x k matrix representing a dictionary of k d-dim. atoms.
%
% outputs:
% DP ..... predicted dictioanry
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$

function DP=mdlsDicPredict(D,method)
if ~exist('method','var')
  method=1;
end

[d,k]=size(D);
w=sqrt(d);
DP=zeros(d,k);

if method==1
  for k1=1:k
    a=D(:,k1);
    ap=zeros(1,d);
    p=0;
    for d1=1:w
      if mod(d1,2) % odd, predict l-r
        for d2=1:w
          idx=(d1-1)*w+d2;
          ap(idx)=a(idx)-p;          
          p=a(idx);
        end
      else % predict r-l
        for d2=w:-1:1
          idx=(d1-1)*w+d2;
          ap(idx)=a(idx)-p;          
          p=a(idx);
        end
      end
    end
    DP(:,k1)=ap;
  end
elseif method==2
  for k1=1:k
    DP(:,k1)=D(:,k1)-sum(D(:,k1))/d;
  end
end