%
% ordering
% creates symbolinc links of the images according to their order in the list
%
% src_dir ..... source folder
% dest_dir .... destination folder
% scheme ...... a name for the classification scheme
% images ..... list of images Must be a Nx1 CELL) where the
%               first column are the filenames without path or extension,
% ext ......... extension of the object files
% aliases ...... cell of Nx1 aliases, one for each image
%
function fname=diskOrder(src_dir, dest_dir, scheme, images, ext, aliases)
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
  %
  % create the links
  %
  printf('Creating links\n');
  for i=1:N
    fname = [ images{i} '.' ext ];
    src_file = [src_dir '/' fname ];
    dest_file = sprintf('%s/%s.%s',dest_dir,aliases{i},ext);
    cmd = sprintf('ln -s %s %s',src_file,dest_file);
    disp(cmd);
    [out,stat]=system(cmd);
  end
