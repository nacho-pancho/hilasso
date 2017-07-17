%%
%% Data caching utility
%% 
%% Given a directory with many object files on it (for example images)
%% the cache controls if there is already computed data for a given object
%% if so, it returns the computed data.
%% 
%% This function loads a cache file
%%
%% inputs:
%% 
%% object_dir ....... Folder where the cached objects live.
%% object_ext ....... Extension of the object files in the folder.
%% cache_file ....... Name of the cache file (without path and extension!)
%%
function cached_data=cache_load(object_dir,object_ext, cache_file)
  if nargin < 3
    error("Must specify at least 3 parameters: dir, extension and cache file name");
  end  
  object_ext=lower(object_ext);

  extindex=rindex(cache_file,'.');
  if !extindex 
    cache_file = [cache_file '.mat'];
  end
  cache_file = [object_dir '/' cache_file];

  if (fexists(cache_file))
    load(cache_file);
  else
    object_list=filenames(object_dir,object_ext);
    N=length(object_list);
    cached_data=cell(N,3); % id, mtime, data
    for i=1:N
      cached_data(i,1)=stripext(object_list(1,i).name);
      cached_data(i,2)=false;
      cached_data(i,3)=[];
    end
    save("-mat",cache_file,"cached_data");
  end
endfunction