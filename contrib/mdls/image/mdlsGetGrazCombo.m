%
% Get labeled data from data annotated database for the list of images requested
% Call with no arguments to get default parameters
% Call only with list_of_images to use default parameters.
%
% INPUT
% list_of_images .... One-dimensional cell where each element contains
%                     the name of an image to be processed (WITH
%                     extension and possibly a RELATIVE path).
%                     You can obtain such a list using mdlsGetFileList.
% params ............ Function parameters.
%
% OUTPUT
% data ........... Normally the patches and labels in a struct. See
%                  mdlsGetData  for a description.
%                  If function was called with no arguments, this is a
%                  struct with the default parameters of the function.
%
function data = mdlsGetGrazCombo(list_of_images, params)
    if nargin < 2
        params = mdlsGetDataAligned();
        params.patch_width = 12;
        params.overlap   = 11;
        params.img_dir   = 'data/inria_graz02/bikes/';
        params.do_align  = false;
        params.do_debug  = false;
        params.remove_dc = true;
        params.save      = true;
        params.do_multiscale = false;
        params.do_sift   = true;
        params.do_patches = true;
        params.mask_window_width = params.patch_width;
        if nargin == 0
            data = params;
            return;
        end
    elseif nargin > 2
        error('Function takes at most 2 arguments.');
    end
    cache_dir = regexprep(params.img_dir,'data','cache');
    system(['mkdir -p ' cache_dir]);
    %
    % weights used for the window analysis. 
    % These are applied to everything: the extracted patches themselves
    % and the mask.
    %
    weights = mdlsCreateGaussianMask(params.patch_width,...
                                     params.mask_window_width);
    weights = weights(:);
    W = sparse(diag(weights)); % weighting matrix
    %
    % we call mdlsGetData or mdlsGetDataAligned depending on params.do_align
    %
    i = 1;
    NI = length(list_of_images);
    NI = length(list_of_images);
    data = struct();
    data.imgSize = [];
    data.gridSize = [];
    data.sift = [];
    data.X    = [];
    data.labels = [];
    data.numPatches = [];
    data.overlap = params.overlap;
    for i=1:NI
        imgpath = list_of_images{i};
        if ~strfind(imgpath,'.image')
            continue;
        end
        matfile = sprintf('%s/%s-combo-w%02d-o%02d-align%d-mscale%d.mat',...
                          cache_dir, ...
                          imgpath(1:end-4),...
                          params.patch_width,...
                          params.overlap,...
                          params.do_align,...
                          params.do_multiscale);
        if exist(matfile,'file')
            fprintf('Loading from %s ...\n', matfile);
            %
            % precomputed!
            %
            load(matfile);
        else
            fprintf('Processing %s ...\n', imgpath);
            if params.do_multiscale
                li = mdlsCreateMultiscaleImages(params.img_dir, ...
                                          list_of_images(i));
            else
                li = list_of_images(i);
            end
            %
            % create multiscale images
            %
            %
            % PATCHES
            %
            if params.do_patches
                if params.do_align
                    thisdata = mdlsGetDataAligned(li,params);
                else
                    thisdata = mdlsGetData(li,params);
                end
                %thisdata.X = single(thisdata.X);
            end
            %
            % SIFT
            %
            if params.do_sift
                patch_size = 32;
                nrml_thres = 1;
                do_padding = 1;
                transposed = true; % we need them in reverse order
                thisdata.numPatches = [];               
                thisdata.sift       = [];
                thisdata.x          = [];
                thisdata.y          = [];
                thisdata.imgSize    = [];
                thisdata.gridSize   = [];
                for ii=1:length(li)
                    % to align the grids, we need to add an extra padding
                    % to the image
                    extra_padding = 32-params.patch_width;
                    I = imread([params.img_dir '/' li{ii}]);
                    I = padarray(I,[extra_padding extra_padding],'replicate','pre');
                    sift = mdlsSift(params.img_dir,...
                                    I,...
                                    grid_spacing,...
                                    patch_size,...
                                    nrml_thres, ...
                                    do_padding,...
                                    transposed);

                    thisdata.numPatches = [thisdata.numPatches size(sift.feaArr,2)];
                    thisdata.imgSize = [thisdata.imgSize ...
                                        [sift.height-extra_padding ; ...
                                        sift.width-extra_padding ]];
                    thisdata.gridSize = [thisdata.gridSize ...
                                        [sift.gridHeight ; sift.gridWidth ]];
                    thisdata.sift = [thisdata.sift sift.feaArr];
                    thisdata.x = [thisdata.x sift.x];
                    thisdata.y = [thisdata.y sift.y];
                end
            end
            %
            % foreground-background binary mask
            %
            maskfname = regexprep(imgpath,'image','mask.all');
            if params.do_multiscale
                li = mdlsCreateMultiscaleImages(...
                    params.img_dir, ...
                    {maskfname});
            else
                li = {maskfname};
            end
            par2=params;
            par2.remove_dc = 0;
            par2.normalize = 0;
            mask = mdlsGetData(li, par2);
            %
            % now we need to assign a label single value to each
            % patch based on the corresponding 1's and 0's in the
            % binary mask
            %
            % we do this by doing a weighted average of the patch
            % using a Gaussian kernel
            %
            thislabel = single(weights' * mask.X);
            save(matfile, 'mask', 'thisdata', 'thislabel');
        end % main loop
        %
        % collate results
        %
        data.imgSize = [data.imgSize thisdata.imgSize];        
        data.gridSize = [data.gridSize thisdata.gridSize];        
        if params.do_patches
            data.X      = [data.X thisdata.X]; clear thisdata.X;
        end
        data.labels = [data.labels thislabel];
        if params.do_sift
            data.sift   = [data.sift thisdata.sift];
            data.x   = [data.x thisdata.x];
            data.y   = [data.y thisdata.y];
        end
        data.numPatches = [data.numPatches thisdata.numPatches];
        if 0 % params.do_debug
            mask2 = double(repmat(labels,params.M,1));
            M = mdlsReconstruct(mask2, mask.imgSize(1,1), mask.imgSize(2,1), ...
                                params.overlap);
            mdlsFigure('Mask'); imagesc(M);            
            clear mask mask2 M;
        end % debug
    end % main loop
end % function