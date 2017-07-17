function D = mdlsGenDCTDictionary(n,K,scramble_dc)
% 
% D = GenerateDCTDictionary(n,K,varargin)
%
% Name: GenerateDCTDictionary
%
% Category:
%
% Description: Create an initial dictionary wit a DCT base.
%
% Input:
% n ........ dimension of the atoms in the dictionary.
% K ........ number of atoms in the dictionary.
%
% Output:
% D ........ the dictionary as a (n x K) matrix.
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%

% % Process the <varargin> arguments
if nargin < 2
    error(['Must specify at leasth dimension and number of atoms. See ' ...
           'help.']);
end
if ~exist('scramble_dc','var')
    scramble_dc = false;
end

% Compute the dictionary
Pn = ceil(sqrt(K));

% Block size is the square root of the dimension
bs = ceil(sqrt(n));
if ~isequal(bs,round(bs))
   fprintf(1,'[ERROR] The dictionary dimension must be a square number.\n')
   D = [];
   return
 end

D = zeros(bs,Pn);
for k = 0:1:Pn-1,
  V = cos((0:1:bs-1)'*k*pi/Pn);
  if k > 0
    V = V - mean(V); 
  end
  D(:,k+1) = V/norm(V);
end
D = kron(D,D);

%
% trim to requested size: discard higher frequency atoms
%
if K < Pn^2
    igv=repmat(1:Pn,1,Pn);
    jgv=repmat(1:Pn,1,Pn);
    jgv=sort(jgv);
    % sort by magnitude and then by angle
    magg = igv.^2 + jgv.^2;
    angg = pi + atan2(jgv,igv);
    sortkey = 1000*2*pi*magg + angg; 
    [gs,gsi] = sort(sortkey);
    idx = gsi(1:K);   
    D = D(:,sort(idx));
end
if scramble_dc
    D(:,1) = randn(size(D,1),1);
    D(:,1) = D(:,1) - mean(D(:,1));
    D(:,1) = D(:,1)./sqrt(sum(D(:,1).^2));
end
