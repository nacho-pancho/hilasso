%
% For each patch in X (where each column is an urolled patch)
% this function returns the 'main direction' of the patch,
% for posterior alignement purposes.
%
function [dir_map,varargout] = mdlsFindMainDirection(X,angles,method,fourier_mask)
if ~exist('angles','var')
    angles = 0:pi/10:pi-pi/10;
end
if ~exist('method','var')
    method = 'fft';
end

%
% create a Gaussian mask for the patch
%
[M,N] = size(X);
w = sqrt(M);
W = mdlsCreateGaussianMask(w,w/4);
%
% apply mask to X
%
WD = spdiags(W(:),0,w^2,w^2);
X = WD*X;
dir_map = zeros(1,N);
if isequal(method,'fft')
    %
    % slow part: find real part of FFT of each patch
    % 
    %wf = 2^ceil(log2(min(8*w,64)));
    wf = w;
    dw=wf-w;
    %if ~exist('fourier_mask','var')
        fourier_mask = mdlsCreateFourierMasks(wf,angles);
        %end
    markers=['|','/','-','\'];
    for j=1:N
        fprintf( '%c', markers(1+mod(floor((j-1)/100),4)) );
        %x = padarray(reshape(X(:,j),w,w),[dw dw],'replicate','post');
        x = reshape(X(:,j),w,w);
        Fx = abs(fft2(x));
        [ma,mai] = max(fourier_mask'*Fx(:));
        dir_map(j) = mai;
        fprintf('\b');
    end
else
    markers=['|','/','-','\'];
    h = fspecial('laplacian');
    w1= 1:w;
    for j=1:N
        fprintf( '%c', markers(1+mod(floor((j-1)/100),4)) );
        x = reshape(X(:,j),w,w);
        x = imfilter(x,h);
        % center of mass
        my = sum((1:w)*x);
        mx = sum((1:w)*x');
        % moments of inertia
        my2 = sum(((1:w)-my).^2*x);
        mx2 = sum(((1:w)-mx).^2*x');
        % angle: atan
        dir_map(j) = 1 + quantiz(mod(atan2(my2,mx2),pi),angles);
        fprintf('\b');
    end
end
if nargout == 2
    varargout{1} = angles;
end
