%%
%% Data caching utility
%% 
%% Given a directory with many object files on it (for example images)
%% the cache controls if there is already computed data for a given object
%% if so, it returns the computed data.
%% 
%% This function shows the status of each entry in the cache:
%% A -1 means the file has changed since its data was cached, 1 means it is up to date
%% and 0 means it has never been computed into the cache.
%% A -2 means there was an unexpected error.
%%
%% inputs:
%% 
%% object_dir ....... Folder where the cached objects live.
%% object_ext ....... Extension of the object files in the folder.
%% cache_file ....... Name of the cache file (without path and extension!)
%%
function data=cache_status(object_dir, object_ext, cached_data)
  if nargin < 3
    error("Must specify at least 3 parameters: dir, extension and cache file name");
  end  
  object_ext=lower(object_ext);

  N=size(cached_data,1);
  data=zeros(N,1);
  for i=1:N
    cache_mtime=cached_data{i,2};
    fname=[object_dir '/' cached_data{i,1} '.' object_ext];
    fstat=stat(fname);
    if (isempty(fstat))
      warning(sprintf("file does not exist: %s",fname));
      data(i,1)=-2;
      continue
    end
    file_mtime=fstat.mtime;
    if fstat.mtime == 0 % never computed
      data(i,1)=0; 
    elseif fstat.mtime > cache_mtime
      data(i,1)=-1;
    else
      data(i,1)=1;
    end
  end
endfunction