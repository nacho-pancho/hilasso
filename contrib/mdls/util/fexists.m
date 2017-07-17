function r=fexists(fname)
  r = !isempty(stat(fname));
endfunction
