%
%
% Input:
% 
% image_dir ...... directory where the images are stored.
% image_list ..... which images are to be processed. Defaults to all PGM
%                  files in image_dir.
% minsz .......... minimum size of scaled images. Defaults to 32
% ratio .......... One dimensional downscaling ratio: each i-th image generated
%                  will be h_ixw_i where h_i =h/ratio^i and
%                  w_i=w/ratio^i, i=1,2,...,i_max where i_max is
%                  arg max_i { min{h_i,w_i} >= minsz }.
%                  Defaults to 2.
%
% Output: nothing as a function, but for each image name.ext, i_max files will 
% be generated with the name being name-{w_i}x{h_i}.ext.
%
% None
% 
function [ms_image_list] = mdlsCreateMultiscaleImages(image_dir, ...
                                                      image_list, ...
                                                      minsz)
    if nargin < 2
        image_list = mdlsGetFileList(image_dir,'pgm');
    end
    if ~exist('minsz','var')
        minsz = 64;
    end        
    if ~exist('ratio','var')
        ratio = 2;
    end

    ratio = 1/ratio;
    ms_image_list = {};
    for i=1:length(image_list)
        ms_image_list{end+1} = image_list{i};
        image_name = [image_dir '/' image_list{i} ];
        fprintf('Processing %5d/%5d',i,length(image_list));
        [pa,na,ex] = fileparts(image_name); 
        I=imread(image_name);
        if size(I,3) > 1
            I=rgb2gray(I);
        end
        [m,n]=size(I);
        s = 1;
        while ratio*min(m,n) > minsz % this is the smallest size we admit
            m = ceil(m*ratio); n = ceil(n*ratio);
            imfile = sprintf('%s-%04dx%04d%s',na,n,m,ex);
            impath = [image_dir '/' imfile];
            if ~exist(impath,'file')
                I = imresize(I, ratio);
                s = ratio*s;
                imwrite(I,impath);
            end
            ms_image_list{end+1} = imfile;
        end
        fprintf('\r');
    end
    fprintf('\n');
end