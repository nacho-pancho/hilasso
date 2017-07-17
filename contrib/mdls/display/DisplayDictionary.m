function varargout = DisplayDictionary(D,varargin)
%
% I = DisplayDictionary(D,margin,mcolor,mag)
%
% Name: DisplayDictionary
%
% Category:
%
% Description: displays a dictionary as a tiled image where each tile is an
% atom assumes square patches, so sqrt(d) must be an integer, where d =
% size(D,1); if there is not output variable the dictionary is shown as an
% image, if there is an output variable the dictionary is returned in the
% first variable.
%
% Input:
% D ........ dictionary to display
% margin ... optional margin to leave between patches. Defaults to 1
% mcolor ... margin color (RGB tern, [0-255]^3). Defaults to [0,0,255]
% (blue)
%
% Output:
% varargout .. an image showing the patches in their original (square) shape,
% all laid out in a grid. An optional magnification can be specified for
% visualization purposes, as patches are usually very small in terms of
% screen resolution.
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%

[d,k] = size(D);
MD = max(max(D));
mD = min(min(D));
sD = 255/(MD-mD);

%
% Process the <varargin> arguments
%
mag = 1;
margin = 5;
mcolor = [0,0,0];
atomUsage = [];
atomChange = [];
plain = 0;
forcedHeight = 0;
for arg = 1:2:length(varargin)
  if isequal(varargin{arg}, 'MarginSize')
    margin = varargin{arg + 1};
  end
  if isequal(varargin{arg}, 'MarginColor')
    mcolor = varargin{arg + 1};
  end
  if isequal(varargin{arg}, 'Magnification')
    mag = varargin{arg + 1};
  end
  if isequal(varargin{arg}, 'AtomUsage')
    atomUsage = varargin{arg + 1};
  end
  if isequal(varargin{arg}, 'AtomChange')
    atomChange = varargin{arg + 1};
  end
  if isequal(varargin{arg}, 'Plain')
    plain = varargin{arg + 1};
  end
  if isequal(varargin{arg}, 'Height')
    forcedHeight = varargin{arg + 1};
  end
end

%
% Coherence
%
if plain == 0
  G = abs(D'*D);
  G = G - diag(diag(G));
  coherence = max(G); % takes only the maximum coherence of each atom w.r.t. the others
else
  coherence = [];
end

% cumulativeCoherence = sum(G)/k; % adds up the coherence of each atom w.r.t. all the others

%
% width of displayed patch, optionally augmented to mag times
%
w = sqrt(d);
% trick to create magnification matrix
magmat = ones(w,1)*floor(1:(1/mag):(w+1)-(1/mag));
cmp = (1:w)'*ones(1,size(magmat,2));
magmat = double(magmat == cmp);
clear cmp;
wm = w*mag;
%
% image tile size, including margin
%
wi = wm + margin;
%
% width and height of tile grid
%
if forcedHeight == 0
  mg = floor(sqrt(k));
else
  mg = forcedHeight;
end
ng = ceil(k/mg);
%
% width and height of output image
%
mi = mg*wi + margin; ni = ng*wi + margin;
I(:,:,1) = mcolor(1)*ones(mi,ni,'uint8');
I(:,:,2) = mcolor(2)*ones(mi,ni,'uint8');
I(:,:,3) = mcolor(3)*ones(mi,ni,'uint8');
i = 0;

for ig = 1:mg
  ii = (ig-1)*wi + 1;
  for jg = 1:ng
    i = i + 1;
    if (i > k)
      break;
    end
    ji = (jg-1)*wi + 1;
    % shape into patch
    patch = reshape(D(:,i),w,w);
    % adjust to 0-255
    patch = round(sD*(patch-mD));
    % if magnification is asked for, enlarge patch:
    if mag > 1
      patch = magmat' * patch * magmat;
    end
    %a=rowCumulativeCoherence(i);
    if plain == 0
      a = coherence(i);    
    else
      a = 0.5;
    end
    if ~isempty(atomUsage)
      tt = 2;
      dofon = 255*atomUsage(i)*ones( (size(patch)+2*tt) );
      I(ii+margin-tt:ii+margin+wm-1+tt,ji+margin-tt:ji+margin+wm-1+tt,:) = cat(3,dofon,dofon,dofon);    
    end
    I(ii+margin:ii+margin+wm-1,ji+margin:ji+margin+wm-1,:) = cat(3,a*patch,0.5*patch,(1-a)*patch);
    if ~isempty(atomChange)
      a = atomChange(i);
      dofon = 255*ones(2,2);
      I(ii+margin-1:ii+margin,ji+margin-1:ji+margin,:) = cat(3,0*dofon,a*dofon,0*dofon);    
    end
  end
  if (i > k)
    break;
  end
end

if isequal(nargout,1)
  % Return the displayable dictionary
  varargout{1} = I;
elseif isequal(nargout,0)
  % Show the displayable dictionary
  imagesc(I);
  colormap(gray);
  axis off
  axis image
end

