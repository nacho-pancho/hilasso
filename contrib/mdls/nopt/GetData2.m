function data = GetData2(patchWidth,overlap,listOfImages)
% 
% data = GetData2(patchWidth,overlap,listOfImages)
%
% Name: GetData2
%
% Category: utility function
%
% Description: extract patches from images
%
% Input:
% patchWidth ..... Width of patches.
% overlap ........ Number of pixels of overlap between the patches (\in
%                  [0,patchWidth-1]).
% listOfImages ... Images to use. This is a cell array where each element
%                  indicates the path to the image file.
%
% Output:
% data ........ Struct with information relating the patches and the images 
%                where they come from. Each struct has the following fields:
% X ........... extracted patches from all images put together as column vectors.
%	DC .......... DC of patches as a row vector
% numImages ... number of images processed
% imgPath ..... list (cell) of path to each image
% imgSize ..... 2xI matrix with the rows and columns of each image
%               as a 2x1 vector.
% numPatches .. 1xI, number of patches in each image
% patchWidth .. width of patches.
% overlap ..... overlap used.                   
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%

nImg = length(listOfImages);
period = patchWidth - overlap;
imgSize = zeros(2,nImg);
numPatches = zeros(1,nImg);
X = [];
DC = [];
for k = 1:nImg 
  % Read the image
  im = double(imread(listOfImages{k}))*(1/255);
  [rows,cols] = size(im);
  imgSize(1,k) = rows; 
  imgSize(2,k) = cols;
  XX = zeros(patchWidth*patchWidth,ceil(rows*cols/period/period));
  DCt = zeros(1,ceil(rows*cols/period/period));

  goRow = true; r = 1;
  goCol = true; c = 1;
  ind = 1;
  while goRow 
    rangeR = r:(r + patchWidth - 1);
    while goCol
      rangeC = c:(c + patchWidth - 1);
      patch = im(rangeR,rangeC);
      temp = patch(:);
      DCt(ind) = mean(temp);
      XX(:,ind) = temp - DCt(ind);
      ind = ind + 1;
      c = c + period;
      if (c > cols - patchWidth + 1), goCol = false; end
    end % while
    r = r + period;     
    goCol = true; c = 1;
    if (r > rows - patchWidth + 1), goRow = false; end
  end % while
  ind = ind - 1;
  X = [X XX(:,1:ind)];
  DC = [DC DCt(1:ind)];
  numPatches(k) = ind;
end % for

% struct
li=cell(1,1); li{1,1}=listOfImages;
data = struct(...
  'X',X,...
  'DC',DC,...
  'numImages',nImg,...
  'imgPath',li,...
  'imgSize',imgSize,...
  'numPatches',numPatches,...
  'patchWidth',patchWidth,...
  'overlap',overlap);

end % function
