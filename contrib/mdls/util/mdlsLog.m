%
% function mylog(fhandle,text)
%
% purpose:
%   print a message to a log file and to the screen
%   simultaneously.
%
% input:
% fhandle ... handle to logfile
% text ...... text to print
%
% output:
% none
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function mdlsLog(fhandle,text)
  fprintf(fhandle,text);
  fprintf(text);
