%
% displays a dictionary as a tiled image where each tile is an atom
% assumes square patches, so sqrt(d) must be an integer, where d=size(D,1)
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
% I ........ an image showing the patches in their original
% (square) shape, all laid out in a grid. An optional
% magnification can be specified for visualization purposes, as
% patches are usually very small in terms of screen resolution.
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function I=mdlsDictDisplay(D,varargin)
%
% parse options
%
margin=1;
mcolor=0;
mag=4;
mwidth=-1;
do_orient = false;
one_dim = false;
do_mosaic = false;
if length(varargin)==1 && iscell(varargin)
    varargin=varargin{1};
end
for i=1:2:length(varargin)
    if isequal(lower(varargin{i}),'margin')
        margin=varargin{i+1};
    elseif isequal(lower(varargin{i}),'mcolor') 
        mcolor = varargin{i+1};
    elseif isequal(lower(varargin{i}),'mag') 
        mag = varargin{i+1};
    elseif isequal(lower(varargin{i}),'mwidth') 
        mwidth = varargin{i+1};
    elseif isequal(lower(varargin{i}),'orient') 
        do_orient = varargin{i+1};
    elseif isequal(lower(varargin{i}),'1d') 
        one_dim = varargin{i+1};
    elseif isequal(lower(varargin{i}),'mosaic') 
        do_mosaic = varargin{i+1};
    else
        error(['Unknown option:' varargin{i}]);
    end
end
%
% multi-dictionary case
%
if iscell(D)
    ND = length(D);
    I=cell(1,ND);
    for i=1:ND
        I{i} = mdlsDictDisplay(D{i},varargin);
    end
    if do_mosaic
        mg = ceil(sqrt(ND));
        ng = ceil(ND/mg);
        [m,n]=size(I{1});
        m2 = m*mg + 5*margin*(mg-1);
        n2 = n*ng + 5*margin*(ng-1);
        dm = m+5*margin;
        dn = n+5*margin;
        I2 = zeros(m2,n2);
        for i=1:mg
            for j=1:ng
                m1=dm*(i-1)+1;
                m2=m1+m-1;
                n1=dn*(j-1)+1;
                n2=n1+n-1;
                i2 = (i-1)*ng+j;
                if i2 > ND
                    break;
                end
                I2(m1:m2,n1:n2) = I{i2};
            end
        end
        I=I2;
    end
    return;
end

if one_dim
    I=mdlsDicDisplayOneD(D,varargin);
    return;
end

[d,k]=size(D);
MD=max(max(D));
mD=min(min(D));
sD=255/(MD-mD);
%
% width of displayed patch, optionally augmented to mag times
%
w=sqrt(d);
% trick to create magnification matrix
magmat=ones(w,1)*floor(1:(1/mag):(w+1)-(1/mag));
cmp=(1:w)'*ones(1,size(magmat,2));
magmat=double(magmat==cmp);
clear cmp;
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

if ~do_orient
    for ig=1:mg
        ii=(ig-1)*wi+1;
        for jg=1:ng
            i=i+1;
            if (i > k)
                break;
            end
            ji=(jg-1)*wi+1;
            % shape into patch
            atom=D(:,i);
            patch=reshape(atom,w,w);
            % adjust to 0-255
            patch=round(sD*(patch-mD));
            % if magnification is asked for, enlarge patch:
            if mag > 1
                patch=magmat' * patch * magmat;
            end
            I(ii+margin:ii+margin+wm-1,ji+margin:ji+margin+wm-1)=patch;
        end
        if (i > k)
            break;
        end
    end
else
    % compute interpolated FFT of each atom
    % the FFT will have the grid size of the tile
    
    Dfft = zeros(mg^2,k);
    % 2D hanning mask
    mask = hann(2*(w+1));
    mask = mask(1:w);
    mask = mask*mask';
    available = 1:k;
    for i=1:k
        atom = reshape(D(:,i),w,w);
        padded_atom = zeros(2*mg,2*mg);
        padded_atom(1:w,1:w) = atom.*mask;
        atom_fft = abs(fft2(padded_atom));
        atom_fft = atom_fft(1:mg,1:mg);
        Dfft(:,i) = reshape(atom_fft,mg^2,1);
    end
    %
    % fill out from the center outwards
    %
    igv=repmat(1:mg,1,ng);
    jgv=repmat(1:ng,1,mg);
    jgv=sort(jgv);
    % sort by magnitude and then by angle
    magg = igv.^2 + jgv.^2;
    angg = pi + atan2(jgv,igv);
    sortkey = 1000*2*pi*magg + angg; 
    [gs,gsi] = sort(sortkey);
    igv = igv(gsi);
    jgv = jgv(gsi);
    for gi=1:length(igv)
        ig = igv(gi);
        jg = jgv(gi);
        ii=(ig-1)*wi+1;
        ji=(jg-1)*wi+1;
        navail = length(available);
        if navail == 0
            break;
        end
        % choose patch from remaining unallocated patches 
        % which has more energy at the corresponding point in the FFT
        [m,i] = max(Dfft((jg-1)*ng+ig,:));
        i2 = available(i);
        % delete just chosen from available list
        available = available([1:(i-1) (i+1):navail]);
        Dfft = Dfft(:,[1:(i-1) (i+1):navail]);
        atom=D(:,i2);
        patch=reshape(atom,w,w);
        % adjust to 0-255
        patch=round(sD*(patch-mD));
        % if magnification is asked for, enlarge patch:
        if mag > 1
            patch=magmat' * patch * magmat;
        end
        I(ii+margin:ii+margin+wm-1,ji+margin:ji+margin+wm-1)=patch;
    end
end

%
% show image
%
if nargout == 0
    imagesc(I); 
    colormap(gray);
end


