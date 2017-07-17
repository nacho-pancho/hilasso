% function mdlsCreateDictImages(dir,marginSize,marginColor,magnification)
%
% Category: utility
%
% Purpose: generate PNG mosaic images of dictionaries found in a directory.
%
%
function mdlsCreateDictImages(dir,varargin)
    close all;
    fileList = mdlsGetFileList(dir,'mat');
    for i=1:length(fileList)
        fname = fileList{i};
        clear D;
        load([dir '/' fname ]);
        if ~exist('D','var')
            continue;
        end
        [path name ext] = fileparts(fname);
        if iscell(D)
            Ds = D; clear D;
            for i=1:length(Ds)
                D=Ds{i};
                fprintf('fname=%s\tclass=%d\trho=%f\n',fname,i,norm(D'*D));
                I=mdlsDictDisplay(D,varargin{:});
                imwrite(I*(255/max(I(:))),sprintf('%s/%s-%d.png',dir,name,i));
            end
        else
            fprintf('fname=%s\trho=%f\n',fname,norm(D'*D));
            I=mdlsDictDisplay(D,varargin{:},'orient',1);
            imwrite(I*(255/max(I(:))),[ dir '/' name '.png']);
        end
    end
end

