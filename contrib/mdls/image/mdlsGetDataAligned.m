function data = mdlsGetDataAligned(list_of_images,params)
% 
% data = mdlsGetData(patch_width,overlap,list_of_images)
%
% Category: utility function
%
% Description: extract patches from many images and concatenate them.
% Call with no arguments to get default parameters
% Call only with list_of_images to use default parameters.
%
% Input:
% list_of_images .... One-dimensional cell where each element contains
%                     the name of an image to be processed (WITH
%                     extension and possibly a RELATIVE path).
%                     You can obtain such a list using mdlsGetFileList.
% params ............ Function parameters. Struct with the following members:
%    patch_width ...... Width of patches.
%    overlap .......... Number of pixels of overlap between the patches (\in
%                       [0,patch_width-1]).
%                       indicates the path to the image file.
%    img_dir .......... Folder where images reside. Defaults to
%    data/images
%
% Output:
%   data ........ Struct with information relating the patches and the images 
%                 where they come from. Each struct has the following fields:
%     X ........... extracted patches from all images put together as column vectors.
%     DC .......... DC of patches as a row vector
%     numImages ... number of images processed
%     imgPath ..... list (cell) of path to each image
%     imgSize ..... 2xI matrix with the rows and columns of each image
%               as a 2x1 vector.
%     numPatches .. 1xI, number of patches in each image
%     patch_width .. width of patches.
%     overlap ..... overlap used.                   
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%
if nargin < 2
    params = struct();
    params.patch_width = 8;
    params.overlap = 7;
    params.img_dir = 'data/images/';
    params.do_debug = false;
    params.remove_dc = true;
    params.force_bw = true;
    params.normalize = true;
    params.method = 'fft';
    params.aprec = 45;
    if nargin == 0
        data = params;
        return;
    end
elseif nargin > 2
    error('Function takes at most 2 arguments.');
end

nImg = length(list_of_images);
imgSize = zeros(2,nImg);
gridSize = zeros(2,nImg);
paddedImgSize = zeros(2,nImg);
numPatches = zeros(1,nImg);
ov = params.overlap;
X=[];
DC=[];
angMap = [];
% shorter names 
w = params.patch_width;
gridSpacing = w-ov;
% size of rotated patches
wr = 2*ceil(w/2*sqrt(2));
dw = wr-w;
dw2 = dw/2;
%
% These indexes are used later for cropping rotated patches
% back to original size
%
tmp = padarray(ones(w,w),[dw2 dw2],0,'both');
cropidx = find(tmp);
clear tmp;

aprec = params.aprec;
angles = 0:(pi/aprec):(pi-pi/aprec);
wf = 2^ceil(log2(min([8*wr 64])));
%
% Fourier mask used to find rotations -- see below mdlsFindMainDirection
%
fourier_mask = sparse(mdlsCreateFourierMasks(wf,angles));
if params.do_debug
    mdlsFigure('Fourier masks');
    mdlsDictDisplay(fourier_mask,'orient',0);
end
%
% Rotation matrices
%
RR = mdlsCreateBatchRotation(wr,180/pi*angles,'bilinear');
%
% main loop over images
%
markers=['|','/','-','\'];
for k = 1:nImg 
    fprintf( '%c', markers(1+mod(k-1,4)) );
    %
    % Read the image  
    %
    im = imread([params.img_dir '/' list_of_images{k}]);
    if (size(im,3) > 1) &&  params.force_bw        
        im = single(rgb2gray(im));
    else
        im = single(im);
    end
    if params.normalize
        im = im*(1.0/255.0);
    end

    [rows,cols] = size(im);
    step = w - ov;
    npr = floor((rows-w)/step)+1;
    npc = floor((cols-w)/step)+1;
    rows2 = (npr-1)*step + w;
    cols2 = (npc-1)*step + w;
    NPt = npc*npr;
    gridSize(:,k) = [npr;npc];
    %
    % to get the rotated patches, we extract bigger patches so that for any
    % rotation we can fill the original ones with real pixels.
    % To get a 1-to-1 correspondence between patches of both sizes, we
    % need to draw the patches from a large image:
    %
    Ir = single(padarray(im,[dw2 dw2],'replicate','both'));
    clear im;
    rowsr = rows + dw;
    colsr = cols + dw;
    %
    % pad if needed so that all patches are well defined
    %
    ovr = ov + dw;
    stepr = wr - ovr;
    nprr = floor((rowsr-wr)/stepr)+1;
    npcr = floor((colsr-wr)/stepr)+1;
    rowsr2 = (nprr-1)*stepr + wr;
    colsr2 = (npcr-1)*stepr + wr;
    Ir = single(padarray(Ir,[rowsr2-rowsr colsr2-colsr],'replicate','post'));
    [Xr,DCr] = mdlsDeconstruct(Ir, wr, ovr);
    clear Ir;
    %
    % number of patches
    %
    NPt = npcr*nprr;
    Xt = zeros(w^2,NPt);
    DCt = zeros(1,NPt);
    numPatches(k) = NPt;
    %
    % find angle of each patch
    %    
    ang_map = mdlsFindMainDirection(Xr,angles,'fft',fourier_mask);
    % add some spatial regularization
    ang_map2 = reshape(ang_map,npcr,nprr);
    ang_map = medfilt2(ang_map2,[5,5],'symmetric');
    ang_map = ang_map(:);
    if params.do_debug
        mdlsFigure('Angle histogram');
        hist(ang_map,36);
        mdlsFigure('Angle');
        imagesc(reshape(ang_map,npcr,nprr)');
    end
    %
    % rotate each patch in Xr
    %
    Xt = mdlsApplyBatchRotation(Xr,RR,ang_map);
    %
    % crop to original size
    %
    Xt = Xt(cropidx,:);
    Xr = Xr(cropidx,:);
    DCt = mean(Xt);
    Xt  = Xt - repmat(DCt,size(Xt,1),1);
    %
    % some heavy debugging
    %
    if params.do_debug
        % 
        % show some patches and the corresponding rotated ones
        %
        mdlsFigure('Alignement test');
        aux = randperm(size(Xr,2));        
        subplot(121); mdlsDictDisplay(Xr(:,aux(1:256)),'orient',0);
        colormap gray; axis off
        title('Original');
        subplot(122); mdlsDictDisplay(Xt(:,aux(1:256)),'orient',0);
        colormap gray; axis off
        title('Aligned');
    end
    %
    % create a rotated version
    %
    X  = [X Xt ];
    DC = [ DC DCt ];
    angMap = [ angMap ang_map' ];
    imgSize(:,k) = [rows; cols];
    paddedImgSize(:,k) = [rows2; cols2];
end % for


% struct
li=cell(1,1); li{1,1}=list_of_images;
data = struct(...
    'X',X,...
    'DC',DC,...
    'angle',angMap,...
    'numImages',nImg,...
    'imgPath',list_of_images,...
    'imgSize',imgSize,...
    'paddedImgSize',paddedImgSize,...
    'numPatches',numPatches,...
    'patchWidth',params.patch_width,...
    'overlap',params.overlap,...
    'gridSize',gridSize,...
    'gridSpacing',gridSpacing);

end % function
