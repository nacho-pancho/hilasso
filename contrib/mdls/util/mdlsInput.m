%
% function v=mdlsInput(caption,def,mode)
%
% input:
% caption ..... text to be shown when asking for parameters
% def ......... default value
% mode ........ like Matlab input, if anything but an 'n' is 
%               specified, the input is assumed to be text.
%
% outputs:
% v ........... value specified or def if enter is pressed
%
%
function v=mdlsInput(caption,def,mode)
  if nargin < 3 
    mode='n';
  end
  if mode=='n'
    v=input(sprintf('%s [%f]:',caption,def));
  else
    v=input(sprintf('%s [%s]:',caption,def),'s');
  end
  if isempty(v)
    v=def;
  end
  
