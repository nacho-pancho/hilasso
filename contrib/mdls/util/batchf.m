%
% Allows the computation of a function on a set of images
% given their type and the folder where they live.
% This script computes the features for each image and stores them
% in a cache file in the same folder so that previously computed images
% do not have to be computed again. 
% It returns a matrix where each row is the result of the function
% for the image (this utility is restricted to functions which return
% row vectors.
%
% This version is for functions (algo_name) which take the full path
% of the image as input. See batch.m for one that takes Matlab matrices as input.
% inputs:
%
% img_dir ....... Folder where the images are located
% img_ext ....... Extension of the images (case insensitive)
% algo_name ..... Name of the function to run on each image. The cache
%                 will be named algo_name.mat
% algo_params ... Optional parameters to pass to the function as row vector
%
function [features,objects]=batchf(img_dir,img_ext,algo_name,algo_params,M)
if nargin < 4
  algo_params=[];
end
%
% load cache
%
cache_data=cache_load(img_dir,img_ext,algo_name);
%
% number of entries in cache
%
N=size(cache_data,1);
objects=cell(N,1);
printf("Processing %d orig images\n",N);    
if N==0 
  warning(sprintf("No images of type %s found in %s.",img_ext, img_dir)); 
  return;
end
if nargin > 4
  if M < N
    N = M;
  end
end

new_data=0;
for i=1:N
  img_name = cache_data{i,1};
  tic;
  entry=cache_read(img_dir, img_ext, cache_data, img_name);
  img_file=[img_dir '/' img_name '.' img_ext];
  if (!isempty(entry))  
    f = entry;
  else
    f = feval(algo_name,img_file,algo_params);
    f = reshape(f,1,prod(size(f))); % only row vectors allowed
    new_data=1;
  end
  time_minute=toc/60; 
  features(i,:) = f;
  objects(i,1) = img_name;
  if i==1, 
    totaltime = round(time_minute*N);
    fprintf("Estimated processing time is %d minutes.\n",totaltime), 
  end
% store in cache
  cache_data=cache_write(img_dir,img_ext,cache_data,img_name,f);  
  if (!mod(i,10) & new_data)
    fprintf('Saving cache.\n');
    cache_save(img_dir,img_ext,cache_data,algo_name);
    new_data = 0;
  end
end
clear X
cache_save(img_dir,img_ext,cache_data,algo_name);
endfunction