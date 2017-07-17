%%
%% Data caching utility
%% 
%% Given a directory with many object files on it (for example images)
%% the cache controls if there is already computed data for a given object
%% if so, it returns the computed data.
%% 
%% This function reads an entry from the cache
%% 
%% inputs:
%% 
%% object_dir ....... Folder where the cached objects live.
%% object_ext ....... Extension of the object files in the folder.
%% cache_file ....... Name of the cache file (without path and extension!)
%% object_id ........ Identifier for a single object in the folder, which
%%                    is nothing but the file name without path or extension.
%%
function data=cache_read(object_dir, object_ext, cached_data, object_id)
  if nargin < 4
    error("Must specify 4 parameters: dir, extension, cache file and object id");
  end  
  object_ext=lower(object_ext);

N=size(cached_data,1);
%
% read
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
  warning("id %s not found in cache: ",object_id);
  data=[];
  return;
end
if cached_data{i,2}>0
  data=cached_data{i,3};
else
  data=[];
end
endfunction