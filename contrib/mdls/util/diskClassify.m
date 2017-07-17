%
% folder-driven class labeling
% creates folders and symbolinc links of the images according
% to their class.
%
% src_dir ..... source folder
% dest_dir .... destination folder
% scheme ...... a name for the classification scheme
% images ..... list of classified images Must be a Nx1 CELL) where the
%               first column are the filenames without path or extension,
% ext ......... Extension of the object files
% labeling .... Vector of labels assigned to each image. Must have same size as images 
%
function fname=diskClassify(src_dir, dest_dir, scheme, images, ext, labeling)
  if nargin != 6
   error("Function accepts exactly 6 arguments");
   return
  end
  [N,dum]=size(images);
  dest_dir = [ dest_dir '/' scheme ];
  [out,stat]=system(['mkdir -p ' dest_dir]);
  if stat != 0
   error(['Could not create ' dest_dir]);
   return
  end
  labels=unique(labeling);
  %
  % create the folders
  %
  L=length(labels);
  for i=1:L
    class_dir = sprintf('%s/%d',dest_dir,labels(i));
    printf('Creating dir %s\n',class_dir);
    [out,stat]=system(['mkdir -p ' class_dir]);
  end
  %
  % create the links
  %
  printf('Creating links\n');
  for i=1:N
    fname = [ images{i} '.' ext ];
    src_file = [src_dir '/' fname ];
    dest_file = sprintf('%s/%d/%s',dest_dir,labeling(i),fname);
    cmd = sprintf('ln -s %s %s',src_file,dest_file);
    disp(cmd);
    [out,stat]=system(cmd);
  end
  %
  % pasting cached data
  %
%  wam_file = [src_dir '/wamFeatures.mat'];
%  if exist(wam_file)
%    %  contains cached_data
%    % which is a cell of Nx3 where first column is image name, 
%    % second is timestamp and 3rd is feature vector
%    load wam_file; 
%  end
%  % split cached_data into one smaller cached_data per label, and save it in dest folder
%  for l=labels
%    li = find(labels==labeling(i));
%    ni = length(li);
%    cached_data2=cell(ni,3);
%    for i=1:ni
%      iname = images{li(i)};
%      tstamp = cached_data{};
%      src_file = [src_dir '/' fname ];
%      dest_file = sprintf('%s/%d/%s',dest_dir,labeling(i),fname);
%      cmd = sprintf('ln -s %s %s',src_file,dest_file);
%      disp(cmd);
%      [out,stat]=system(cmd);
%    end
%   end
endfunction
