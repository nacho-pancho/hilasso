%
% Return date in yyyy-mm-dd:HH:MM:SS format
% function s=mdlsDateString()
%
%
function s=mdlsDateString()
    s=datestr(now,'yyyy-mm-dd-HH.MM.SS');
end