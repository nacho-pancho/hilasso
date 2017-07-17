%
% Generic grid search implementation
% given a function name (or handle) and a cell with arrays of parameter values,
% calls the function with each combination of parameters and stores the
% result in a multi-dimensional array.
%
% INPUT
% fun .......... (string or handle) specifies the function. 
%                The function takes a number of arguments and
%                returns a scalar value. The numbre of argumenst is
%                determined by funargs (see below).
% funargs ...... Array of structs, one per parameter. 
%                Each struct has the following fields:
%   name ......... name of parameter
%   values ....... cell with values to test
%
function G = gengrid(fun,funargs)

if ischar(fun)
  eval(['fun=@%s' fun]);
end

scriptname = sprintf('gengrid_tmp_%d.m', int32(ceil(rand()*1e5)))
NARGS=length(funargs);
argnames = cell(1,NARGS);
argvars = cell(1,NARGS);
argvalues = cell(1,NARGS);
%
% create a Matlab script from the parameters
%
fh = fopen(scriptname,'w');
indent = 0;
callstr = 'fun(';
resstr = 'res(';
resinistr = 'res = zeros(';
for i=1:NARGS
    argnames{i}     = funargs(i).name;
    argvalues{i}     = funargs(i).values;
    idx_name = [argnames{i} '_idx'];
    fprintf(fh,'%sfor %s=1:length(argvalues{par_i})\n',...
            repmat(' ',1,indent), ...
            idx_name);
    indent = indent + 2;
    callstr = [callstr 'argvalues{par_i}{' idx_name '}'];
    resstr  = [resstr  idx_name];
    resinistr  = [resinistr  num2str(length(argvalues{i}))];
    if i == NARGS 
        callstr = [callstr ')'];
        resstr = [resstr ')'];
        resinistr = [resinistr ')'];
        %
        % inner loop: call function
        %
        fprintf(fh,'%s%s = %s\n',repmat(' ',1,indent),resstr,callstr);
    else
        fprintf(fh,'%spar_i = par_i + 1;\n',...
                repmat(' ',1,indent));
        callstr = [callstr ', '];
        resstr = [resstr ', '];
        resinistr = [resinistr ', '];
    end
end
for i=1:NARGS
    indent = indent - 2;
    fprintf(fh,'%send\n',repmat(' ',1,indent));
    fprintf(fh,'%spar_i = par_i - 1;\n',...
            repmat(' ',1,indent));
end
fclose(fh);
%
% call dynamically-created script
%
par_i = 1;
eval(resinistr);
eval(['run ./' scriptname]);
%
% delete it
%
%delete(scriptname);
end