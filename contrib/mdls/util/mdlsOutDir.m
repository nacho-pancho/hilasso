%
% Possibly creates and returns the output directory corresponding  to an experiment for the specificed date (current by default).
% If the date is specified, it is assumed that we want to READ an
% existing set of results, so nothing is created.
%
function outdir=mdlsOutDir(exp_name,when)
    if nargin == 0
        error('Must specify experiment name.');
    end
    if nargin < 2
        % write mode
        when = now();
        ds = datestr(when,'yyyy-mm-dd-HH.MM');
        outdir = [ 'results/' exp_name];
        if ~exist(outdir,'file')
            mkdir('.', outdir);
        end
        outdir = [ outdir '/' ds ];
        if ~exist(outdir,'file')
            mkdir('.', outdir);
        end
    else
        % read mode
        if isequal(when,'latest')
            ds='latest';
        else
            ds = datestr(when,'yyyy-mm-dd-HH.MM');
        end
        outdir = [ 'results/' exp_name '/' ds];
    end    
end