%
% purpose: 
% load an image given its bare name. To be used in MDL-S project.
% requires file mdlsConstants.m to be in search path.
%
% input:
% imgname .... image name with no path or extension, e.g. 'boat'.
%
% output:
% I .......... the image as a full double matrix
%
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function I=mdlsLoadImage(imgname,normalize)
  mdlsConstants;
  I=double(imread([IMAGES '/' imgname '.png']));
  if exist('normalize','var') && normalize
    I=I/255;
  end
end