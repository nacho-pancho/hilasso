function [Xr,ang] = mdlsAlignPatches(X,varargin)
% 
% data = mdlsGetData(patchWidth,varargin (as key-value pairs))
%
% Category: utility function
%
% Description: take a set of patches and align them so that the main direction points east
%
% Input:
% X ..................... unaligned patches.
% angular_prec (key) .... angular precision (opt). Defaults to pi/180
% fourier_over (key) .... oversampling factor for Fourier transform. Defaults to 
% 
% Output:
% Xa ........... aligned patches
% ang .......... angles of the patches
%
% Author: I. Ramirez <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%
ang_prec = pi/180;
do_debug = 0;
fourier_over = 8;
interp_method = 'bicubic';
max_rotation = 0;
for i=1:2:length(varargin)
    if isequal(varargin{i},'max_rotation')
        max_rotation = varargin{i+1};
    elseif isequal(varargin{i},'angular_prec')
        ang_prec = varargin{i+1};
    elseif isequal(varargin{i},'fourier_over')
        fourier_over = varargin{i+1};
    elseif isequal(varargin{i},'do_debug')
        do_debug = varargin{i+1};
    elseif isequal(varargin{i},'interp_method')
        interp_method = varargin{i+1};
    else
        error(['Unknown parameter: ' varargin{i}]);
    end
end
[M,N] = size(X);
% shorter names 
w = sqrt(M);
%
% hanning mask used for rotation
%
x_mask = zeros(w,w);
w2 = w/2;
h = hann(10*w2);
for i = 1:w
    for j = 1:w
        r = sqrt((i-0.5-w2)^2+(j-w2-0.5)^2);
        if r < w2
            x_mask(i,j) = h(round(5*(w2+r)));
        end
    end
end
if do_debug
    mdlsFigure('XMask'); imagesc(x_mask);
end
x_mask = x_mask(:);
%
% mask  used to integrate over rays in the fourier domain
%
if max_rotation == 0
    angles = 0:ang_prec:(pi-ang_prec);
else
    %angles = [0:ang_prec:(max_rotation-ang_prec) (pi-max_rotation):(pi-ang_prec)];
    angles = 0:ang_prec:(max_rotation-ang_prec);
end
wf = 2^ceil(log2(min([fourier_over*w 128])));
fourier_mask = sparse(mdlsCreateFourierMasks(wf,angles));
if do_debug
    mdlsFigure('Fourier masks');
    mdlsDicDisplay(fourier_mask,'orient',0);
end
xpad = zeros(wf,wf);
ang = zeros(1,N);
Xr=zeros(size(X));
for j=1:N
    %
    % find maximum radial integral
    %
    xr = reshape(X(:,j).*x_mask,w,w);
    xpad(1:w,1:w) = xr;
    Fx = abs(fft2(xpad));
    radon = Fx(:)'*fourier_mask;
    [ma,mai] = max(radon);
    %
    %
    % now rotate -ang, where ang is the angle with maximum
    % power. The angle computed above is actually clockwise, 
    % while imrotate uses counter-clockwise, so the '-' is implicit.
    %
    ang(j) = angles(mai);
    xr = imrotate(reshape(X(:,j),w,w),ang(j)*180/pi,interp_method,'crop');
    %dc = mean(xr(:));
    %xr = xr - dc;
    Xr(:,j) = xr(:);
end % for

end % function
