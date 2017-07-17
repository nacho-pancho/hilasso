%%
%% return a detector struct with the given data
%%
function det=buildDetector(det_fun,det_params)
  det=struct('fun',det_fun,'params',det_params);
endfunction