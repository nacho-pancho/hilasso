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
% img_dir ....... Folder where the images are located
% img_ext ....... Extension of the images (case insensitive)
% algo_name ..... Name of the function to run on each image. The cache
%                 will be named algo_name.mat
% algo_params ... Optional parameters to pass to the function as row vector
%
%
function [features,objects]=batch(img_dir,img_ext,algo_name,M)
%
% load cache
%
cache_data=cache_load(img_dir,img_ext,algo_name);
N=size(cache_data,1);
objects=cell(N,1);
%
% number of entries in cache
%
if N==0 
  warning(sprintf("No images of type %s found in %s.",img_ext, img_dir)); 
  return;
end

if nargin > 3
  if M < N
    N = M;
  end
end
printf("Processing %d orig images\n",N);    
if nargout > 1
  imglist=cache_data(:,1);
end
new_data=0;
for i=1:N
  img_name = cache_data{i,1};
  entry=cache_read(img_dir, img_ext, cache_data, img_name);
  img_file=[img_dir '/' img_name '.' img_ext];
  if (!isempty(entry))  
    time_minute=0;
    f = entry;
  else
    % first time: have to compute
    printf("%3d/%3d:%s\n",i,N,img_file);
    X = imread(img_file);                 % read one image at a time
    X = double(X);       
    tic
    f = feval(algo_name,X);
    time_minute=toc/60; 
    f = reshape(f,1,prod(size(f))); % only row vectors allowed
    new_data = 1; % we have new data to store
  end
  objects(i,1) = img_name;
  features(i,:) = f;
  if i==1, 
    totaltime = round(time_minute*N);
    fprintf("Estimated processing time is %d minutes.\n",totaltime);
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
