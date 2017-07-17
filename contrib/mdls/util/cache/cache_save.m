%%
%% Data caching utility
%% 
%% Given a directory with many object files on it (for example images)
%% the cache controls if there is already computed data for a given object
%% if so, it returns the computed data.
%% 
%% With the first three arguments a status vector is returned where -1 means
%% the file has changed since its data was cached, 1 means it is up to date
%% and 0 means it has never been computed into the cache.
%%
%% With four arguments, it reads the cache file and looks for the
%% row containing precomputed data for object_id.
%% If no precomputed data is found, it returns an empty vector.
%%
%% if a fifth argument is passed, the cache entry for object_id 
%% is written with new_data and new_data is returned.
%%
%%
%% inputs:
%% 
%% object_dir ....... Folder where the cached objects live.
%% object_ext ....... Extension of the object files in the folder.
%% cache_file ....... Name of the cache file (without path and extension!)
%% object_id ........ Identifier for a single object in the folder, which
%%                    is nothing but the file name without path or extension.
%% new_data ......... If present, write the cache entry for object_id with this
%%                    data. 
%%
%%
%%
function cache_save(object_dir,object_ext, cached_data, cache_file)
  if nargin < 4
    error("Must specify at least 4 parameters: dir, extension and cache file name");
  end  

  extindex=rindex(cache_file,'.');
  if !extindex 
    cache_file = [cache_file '.mat'];
  end
  cache_file = [object_dir '/' cache_file];

  save("-mat",cache_file,"cached_data");

endfunction