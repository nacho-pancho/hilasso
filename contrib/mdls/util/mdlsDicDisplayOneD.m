%
% displays a dictionary as a tiled image where each tile is an atom
% this version is for 1D signals, so a plot is shown in each tile.
%
% inputs:
% D ........ dictionary to display
% margin ... optional margin to leave between patches. Defaults to 1
% mcolor ... margin color (must be scalar between 0-255). Defaults to 1
%            (white)
% width .... if specified, sets the width (in atoms) of the
%            resulting mosaic.
%
% outputs:
% I ........ an image showing the atoms as 1D signals 
%             all laid out in a grid. An optional
% magnification can be specified for visualization purposes, as
% patches are usually very small in terms of screen resolution.
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function I=mdlsDicDisplayOneD(D,varargin)
%
% parse options
%
margin=1;
mcolor=0;
mag=1;
mwidth=-1;
% DEACTIVATED FOR NOW
I=zeros(10,10);
return;

do_orient = true;
if length(varargin)==1 && iscell(varargin)
    varargin=varargin{1};
end
for i=1:2:length(varargin)
    if isequal(lower(varargin{i}),'1d')
        % ignore this one: may come from DicDisplay
    elseif isequal(lower(varargin{i}),'margin')
        margin=varargin{i+1};
    elseif isequal(lower(varargin{i}),'mcolor') 
        mcolor = varargin{i+1};
    elseif isequal(lower(varargin{i}),'mag') 
        mag = varargin{i+1};
    elseif isequal(lower(varargin{i}),'mwidth') 
        mwidth = varargin{i+1};
    elseif isequal(lower(varargin{i}),'orient') 
        do_orient = varargin{i+1};
    else
        error(['Unknown option:' varargin{i}]);
    end
end

[d,k]=size(D);
MD=max(max(D));
mD=min(min(D));
sD=255/(MD-mD);
%
% width of displayed patch, optionally augmented to mag times
%
w=size(D,1);
% trick to create magnification matrix
wm=w*mag;
%
% image tile size, including margin
%
wi=wm+2*margin; 
%
% width and height of tile grid
%
if mwidth <= 0
    mg=floor(sqrt(k));
else
    mg=mwidth;
end
ng=ceil(k/mg);
% override grid size for oriented mode: square
% may contain more tiles than actual 
if do_orient
    mg=ceil(sqrt(k));
    ng = mg;
end

%
% width and height of output image
%
mi=mg*wi; ni=ng*wi;
I=uint8(mcolor*ones(mi,ni));
i=0;

for ig=1:mg
    ii=(ig-1)*wi+1;
    for jg=1:ng
        i=i+1;
        if (i > k)
            break;
        end
        ji=(jg-1)*wi+1;
        % plot atom in a wxw graph
        atom=D(:,i);
        tile = zeros(wm,wm);
        ymin = min(atom);
        ymax = max(atom);
        yscal = (wm-1)/(ymax-ymin);
        xscal = (w-1)/(wm-1);
        for x=1:wm
            xp = 1+xscal*(x-1);
            x0 = floor(xp);
            x1 = x0 + 1;
            if x1 <= w
                y1 = 1 + yscal*(atom(x1)-ymin);
                y0 = 1 + yscal*(atom(x0)-ymin);
                yp = y0*(x1-xp) + y1*(xp-x0);
            else
                yp = 1 + yscal*(atom(x0)-ymin);
            end
            y0 = floor(yp);
            y1 = y0 + 1;
            if y1 <= wm
                tile(wm-y0+1,x) = tile(wm-y0+1,x) + (y1-yp);
                tile(wm-y1+1,x) = tile(wm-y1+1,x) + (yp-y0);
            else
                tile(y0,x) = tile(y0,x) + 1;
            end
        end
        I(ii+margin:ii+margin+wm-1,ji+margin:ji+margin+wm-1)=tile;
    end
    if (i > k)
        break;
    end
end

%
% show image
%
if nargout == 0
    imagesc(I); 
    colormap(gray);
end


