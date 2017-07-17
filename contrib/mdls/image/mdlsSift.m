function feaSet = mdlsSift(img_dir,...
                           I, ...                           
                           grid_spacing, ...
                           patch_width, ...
                           nrml_threshold, ...
                           do_padding, ...
                           transposed)
    %disp('Extracting SIFT features...');

    siftLens = [];
    if ~exist('do_padding','var')
        do_padding = true;
    end

    if iscell(I)
        NI = length(I);
        feaSet = struct();
        feaSet.feaArr = [];
        feaSet.x = [];
        feaSet.y = [];
        feaSet.width = [];
        feaSet.height= [];
        for i=1:length(I)
            f=mdlsSift(img_dir,...
                       I{i},...
                       grid_spacing,...
                       patch_width,...
                       nrml_threshold, ...
                       do_padding,...
                       transposed);
            feaSet.feaArr = [ feaSet.feaArr  f.feaArr ];
            feaSet.x      = [ feaSet.x ; f.x ];
            feaSet.y      = [ feaSet.y ; f.y ];
            feaSet.width  = [ feaSet.width ; f.width ];
            feaSet.height = [ feaSet.height ; f.height ];
        end
        return;
    end

    if ischar(I)
        I = imread([img_dir '/' I]);
    end

    if ndims(I) == 3,
        I = im2double(rgb2gray(I));
    else
        I = im2double(I);
    end

    [im_h, im_w] = size(I);

    if do_padding
        %
        % pad image
        %
        step = grid_spacing;
        w    = patch_width;
        rows = im_h;
        cols = im_w;
        npr = ceil((rows-w)/step)+1;
        npc = ceil((cols-w)/step)+1;
        rows2 = (npr-1)*step + w;
        cols2 = (npc-1)*step + w;
        I = padarray(I,[rows2-rows cols2-cols],'replicate','post');
    end
    [im_h, im_w] = size(I);
    %
    % make grid sampling SIFT descriptors
    %
    remX = mod(im_w-patch_width,grid_spacing);
    offsetX = floor(remX/2)+1;
    remY = mod(im_h-patch_width,grid_spacing);
    offsetY = floor(remY/2)+1;

    grid_x = offsetX:grid_spacing:im_w-patch_width+1;
    grid_y = offsetY:grid_spacing:im_h-patch_width+1;
    grid_n = length(grid_x);
    grid_m = length(grid_y);
    [gridX,gridY] = meshgrid(grid_x, grid_y);

    % find SIFT descriptors
    siftArr = sp_find_sift_grid(I, gridX, gridY, patch_width, 0.8);
    [siftArr, siftlen] = sp_normalize_sift(siftArr, nrml_threshold);

    siftArr=siftArr'; % put descriptors as columns
    if transposed
        gridX = gridX';
        gridY = gridY';
        idx = 1:size(siftArr,2);
        idx = reshape(idx,grid_m,grid_n); 
        idx = idx'; idx = idx(:);
        siftArr = siftArr(:,idx);
    end
    feaSet.feaArr = siftArr;
    feaSet.x = gridX(:) + patch_width/2 - 0.5;
    feaSet.y = gridY(:) + patch_width/2 - 0.5;
    feaSet.gridWidth = grid_n;
    feaSet.gridHeight = grid_m;
    feaSet.width = im_w;
    feaSet.height = im_h;
end