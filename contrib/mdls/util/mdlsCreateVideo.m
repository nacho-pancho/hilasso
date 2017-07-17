%
% Create a video with the evolution of a dictioanry during
% training. Requires the tools "convert" (ImageMagick suite) 
% and ffmpeg and works only in Unix.
%
function CreateVideo(dir,prefix)
 %
 % convert to JPG
 %
 system(sprintf('cd %s; for i in %s*.png; do convert -quality 100 $i ${i/png/jpg}; done',dir,prefix));
 %
 % encode to MPEG
 %
 system([ sprintf('cd %s',dir) '; ffmpeg -intra -r 25 -i d%05d.jpg dic.avi']);
 %
 % cleanup 
 %
 system(sprintf('cd %s; rm -f d*.jpg',dir));

 
