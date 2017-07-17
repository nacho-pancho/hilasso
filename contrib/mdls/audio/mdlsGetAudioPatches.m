function data = mdlsGetAudioPatches(patchWidth,overlap,list_of_tracks, audio_dir, remove_dc)
% 
% function data = mdlsGetAudioPatches(patchWidth,overlap,list_of_tracks,audio_dir, remove_dc)
%
% Category: utility function
%
% Description: extract patches from audio files
%
 % Input:
% patch_width ..... Width of patches.
% overlap ........ Number of pixels of overlap between the patches (\in
%                  [0,patchWidth-1]).
% list_of_tracks . Tracks to use. This is a cell array where each element
%                  indicates the path to an audio file.
% audio_dir ........ Base directory for images, defaults to data/images.
% remove_dc ...... if true, removes the DC from each patch and stores it in the DC variable
%                  otherwise DC is not removed and data.DC will be empty. Default: true.
% Output:
% data ........ Struct with information relating the patches and the images 
%                where they come from. Each struct has the following fields:
%    X ........... extracted patches from all images put together as column vectors.
%    DC .......... DC of patches as a row vector, if remove_dc = true, otherwise empty.
%    numTracks ... number of images processed
%    filePath ..... list (cell) of path to each track
%    fileSize .... 2xI matrix with the length of the file in samples.
%    numPatches .. 1xI, number of patches in each track
%    patchWidth .. width of tracks
%    overlap ..... overlap used.                   
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%
if ~exist('remove_dc','var')
    remove_dc = true;
end
if ~exist('track_dir','var')
    track_dir = 'data/images/';
end

nImg = length(list_of_tracks);
imgSize = zeros(2,nImg);
paddedTrackSize = zeros(2,nTrack);
numPatches = zeros(1,nTrack);
% shorter names 
w = patchWidth;
ov = overlap;
X=[];
DC=[];
for k = 1:nTrack 
  % Read the image
  [ipath,iname,iext] = fileparts(list_of_tracks{k});
  fprintf('loading %s\n',iname);
  if isempty(ipath)
      ipath = track_dir;
  end
  if isempty(iext)
      iext = '.wav';
  end
  lab = [ipath,iname,'.lab'];
  if ~isempty(labfile)
      labdata = load(labfile);
  end
  [im,fs,bs] = single(wavread([ipath '/' iname iext]));
  [nsamp] = length(im);
  %
  % pad if needed so that all patches are well defined
  %
  % determine needed dimension
  %
  step = w - ov;
  np = ceil((nsamp-w)/step)+1;
  np2 = (np-1)*step + w;

  im2 = [ im, zeros(1,np2-np) ];
  %
  % decompose: use hamming or not?
  %
  Xt = zeros(w,np2);
  for si=1:step:np
      Xt(:,si) = im(si:(si+w-1));
  end
  if remove_dc
      DCt = mean(Xt);
      DC = [ DC DCt ];
      Xt = Xt - repmat(DCt,w,1);
  end
  X = [X Xt ];
  trackSize(:,k) = np2;
end % for


% struct
li=cell(1,1); li{1,1}=list_of_tracks;
data = struct(...
  'X',X,...
  'DC',DC,...
  'numImages',nTrack,...
  'trackPath',li,...
  'trackSize',trackSize,...
  'numPatches',numPatches,...
  'patchWidth',patchWidth,...
  'overlap',overlap);

end % function
