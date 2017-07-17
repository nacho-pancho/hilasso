function names=filenames(directory,extension)
% The function lists all filenames of the specified extension from the directory 
% directory  the path to the files to list in names; 
%            if not specified or is '' then the current directory is taken
% extension   can be any file extension with any number of letters;
%        is not case sensitive; 
%        if not specified then all files are listed in the output names
% Typical usage:   names=filenames('C:\My Documents\MATLAB\ImagesDB','jpg');
%                  names=filenamesc('','bmp');

if nargin<1, 
  directory='.'; 
end
if nargin<2
  extension='*';
end

%
% cached?
%
if extension!='*'
  cache_file=[directory '/dirlist-' extension '.mat'];
else
  cache_file=[directory '/dirlist-all.mat'];
end

if size(stat(cache_file),1) > 0
  printf("dirlist found.\n",directory);
  load(cache_file);
  return;
end

A=dir([directory '/']);
numfiles=length(A);
if nargin<2
   for i=3:numfiles
      names(i-2).name=A(i).name;
   end
else
   dotextension=strcat('.',lower(extension));
   let=length(extension);
   counter=0;
   for i=1:numfiles
      nextname=A(i).name;
      if length(nextname)>let
         if lower(nextname(end-let:end))==dotextension,
            counter=counter+1; names(counter).name=nextname;
         end
      end
   end
end

save("-mat",cache_file,"names");
endfunction