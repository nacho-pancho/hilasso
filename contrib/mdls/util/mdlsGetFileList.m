%
% Given a folder, return all the images with the given extension
% as a cell array.
%
% Input: 
% folder ......... the folder where images are stored. This is
%                  non-recursive.
% extension ...... File extension. Defaults to png
%
function flist = mdlsGetFileList(folder,extension,pattern)
    if ~exist('extension','var')
        extension = 'png';
    end
    if ~exist('pattern','var')
        pattern=[];
    end
    fdir=dir(folder);
    flist=cell(1,length(fdir));
    j=1;
    for i=1:length(fdir)
        fname = fdir(i).name;
        dotidx = findstr(fdir(i).name,'.');
        if isempty(dotidx)
            continue;
        end
        fext = fname(dotidx(end)+1:end);
        if isequal(lower(fext),lower(extension))
            if isempty(pattern) || ~isempty(strfind(fname,pattern))
                flist{1,j} = fname;
                j=j+1;
            end
        end
    end
    flist=flist(1:j-1);
end