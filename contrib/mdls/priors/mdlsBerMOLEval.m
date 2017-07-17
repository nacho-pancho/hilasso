%
% purpose:
% evaluate the Bernoulli-Gamma-Laplacian distribution at the given points.
%
% inputs:
% x ........ evaluation point(s)
% theta .... Bernoulli parameter, i.e., p(Z=1)
% k ........ shape of Gamma
% b ........ scale of Gamma
%
% outputs:
% fx ....... the distribution evaluated at x
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function fx=mdlsBerMOLEval(x,theta,k,b)
 fx=theta.*mdlsMOLEval(abs(x),k,b)/2; % divide by two since gamlapPDF is one-sided
 i0=find(x==0);
 if ~isempty(i0)
     fx(i0)=fx(i0)+(1-theta);
 end


