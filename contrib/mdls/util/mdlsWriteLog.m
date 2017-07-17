%
% function WriteLog(S,f,...)
%
% description: writes strings simultaneously to screen and to a text file.
% the file ID is mantained internally and controlled using the special commands
% 'open' and 'close'.
%
% inputs:
% 
%  S ......... printf format string or the reserved strings 'open' or 'close'
%  f ......... file name for 'open', ignored for 'close' or first argument to printf format string
%  (var) ..... further arguments to the format string
% 
% outputs:
%  n ......... number of chars written
%
