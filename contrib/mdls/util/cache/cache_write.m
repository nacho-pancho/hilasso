%%
%% Data caching utility
%% 
%% Given a directory with many object files on it (for example images)
%% the cache controls if there is already computed data for a given object
%% if so, it returns the computed data.
%% 
%% This function writes data into the cache structure
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
function cached_data=cache_write(object_dir,object_ext, cached_data, object_id, new_data)
  if nargin < 5
    error("Must specify at least 3 parameters: dir, extension and cache file name");
  end  
  if nargout != 1
    error("Bad usage: must assign result back to cached_data in order to update it!");
  end
  object_ext=lower(object_ext);
N=size(cached_data,1);
%
% read/write
%
% strip extension from object id if redundant
  extindex = rindex(object_id,'.');
  fext=lower(substr(object_id,extindex+1));
  if strcmp(fext,object_ext)
    object_id = substr(object_id,1,extindex-1);
  end
  found=0;
  for i=1:N
    if (strcmp(cached_data{i,1},object_id))
      found=1;
      break
    end
  end
  if !found
    warning("id %s not found in cache!",object_id);
    data=[];
    return;
  end
  fname=[object_dir '/' object_id '.' object_ext];
  fstat = stat(fname);
  if (length(fstat)==0)
    error(['File ' fname ' not found but is listed in the cache!']);
  end
  cached_data{i,2}=fstat.mtime;
  cached_data{i,3}=new_data;
endfunction