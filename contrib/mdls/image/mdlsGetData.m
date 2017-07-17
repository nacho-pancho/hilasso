function data = mdlsGetData(list_of_images, params)
% 
% function data = mdlsGetData(list_of_images, params)
%
% Category: utility function
%
% Description: extract patches from many images and concatenate them.
%
% Input:
% list_of_images ... Images to use. This is a cell array where each element
%                  indicates the path to the image file.
% params
%     patch_width ..... Width of patches.
%     overlap ........ Number of pixels of overlap between the patches (\in
%                  [0,patchWidth-1]).
%     img_dir ........ Base directory for images, defaults to data/images.
%     remove_dc ...... if true, removes the DC from each patch and stores it in the DC variable
%                  otherwise DC is not removed and data.DC will be
%                  empty. Default: true.
%
% Output:
% data ........ Struct with information relating the patches and the images 
%                where they come from. Each struct has the following fields:
%    X ........... extracted patches from all images put together as column vectors.
%    DC .......... DC of patches as a row vector, if remove_dc = true, otherwise empty.
%    numImages ... number of images processed
%    imgPath ..... list (cell) of path to each image
%    imgSize ..... 2xI matrix with the rows and columns of each image
%                  as a 2x1 vector.
%    numPatches .. 1xI, number of patches in each image
%    patchWidth .. width of patches.
%    overlap ..... overlap used.                   
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
    if nargin == 0
        data = params;
        return;
    end
elseif nargin > 2
    error('Function takes at most 2 arguments.');
end
if ischar(list_of_images)
  list_of_images = {list_of_images};
end

nImg = length(list_of_images);
imgSize = zeros(2,nImg);
paddedImgSize = zeros(2,nImg);
numPatches = zeros(1,nImg);
% shorter names 
w = params.patch_width;
ov = params.overlap;
gridSpacing = w-ov;
X=[];
DC=[];
for k = 1:nImg 
  % Read the image
  [ipath,iname,iext] = fileparts(list_of_images{k});
  %fprintf('loading %s\n',iname);
  if isempty(ipath)
      ipath = params.img_dir;
  end
  if isempty(iext)
      iext = '.png';
  end
  im = imread([ipath '/' iname iext]);
  if (size(im,3) > 1) &&  params.force_bw        
      im = single(rgb2gray(im));
  else
      im = single(im);
  end
  if params.normalize
      im = im*(1.0/255.0);
  end
  [rows,cols] = size(im);
  %
  % pad if needed so that all patches are well defined
  %
  % determine needed dimension
  %
  step = w - ov;
  npr = ceil((rows-w)/step)+1;
  npc = ceil((cols-w)/step)+1;
  rows2 = (npr-1)*step + w;
  cols2 = (npc-1)*step + w;

  im2 = padarray(im,[rows2-rows cols2-cols],'replicate','post');
  %
  % number of patches
  %
  NPt = npc*npr;
  gridSize(:,k) = [npr;npc];
  if params.remove_dc
      [Xt,DCt] = mdlsDeconstruct(im2,w,ov);
  else
      Xt = mdlsDeconstructFast(im2,w,ov);
  end          
  numPatches(k) = NPt;
  X  = [X Xt ];
  if params.remove_dc
      DC = [ DC DCt ];
  end
  imgSize(:,k) = [rows; cols];
  paddedImgSize(:,k) = [rows2; cols2];
end % for


% struct
li=cell(1,1); li{1,1}=list_of_images;
data = struct(...
  'X',X,...
  'DC',DC,...
  'numImages',nImg,...
  'imgPath',li,...
  'gridSize',gridSize,...
  'imgSize',imgSize,...
  'numPatches',numPatches,...
  'patchWidth',params.patch_width,...
  'gridSpacing',gridSpacing,...
  'overlap',params.overlap);

end % function
