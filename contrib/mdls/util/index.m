%
% wrapper for Matlab, since the function is called index in Octave
% 
function idx=index(s,r)
  idx=findstr(s,r);
end