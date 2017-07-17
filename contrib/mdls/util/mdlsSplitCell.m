function [c]=mdlsSplitCell(s,d)
%
% similar to split() but returns a cell of string objects
%
if nargin < 2
 d = ' ';
end

s1=s;
c={};
while ~isempty(s1)
  idx = index(s1,d);
  if isempty(idx)
    c={c{:} s1};
    return
  elseif idx == 1
    continue
  else
    c={c{:} s1(1:idx-1) };
    s1=s1(idx+1:end);
  end
end

