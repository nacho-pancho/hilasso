function [s]=mdlsJoinCell(c,d)
%
% similar to join() but works with a cell of string objects
%
if nargin < 2
 d = ' ';
end

if length(c) == 0
  s='';
  return;
else
  s=c{1};
end
  
for i=2:length(c)
 s=[s d c{i}];
end
