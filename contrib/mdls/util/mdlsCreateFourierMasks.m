% 
% function fmask = mdlsCreateFourierMask()
%
% Category: surrogate function
%
% Description: used by mdlsGetData to align patches prior
%
% Input:
%
% 
%
% Output:
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Id$
%
function fourier_mask = mdlsCreateFourierMasks(w,angles)
    w2 = (w-1)/2;
    fourier_mask=zeros(w^2,length(angles));
    h = hann(10*w2);
    for ai = 1:length(angles)
        ang = angles(ai);
        mask=zeros(w,w);
        for r = 1 : 0.1 : w2-1
            i = r*sin(ang);
            j = r*cos(ang);
            %fprintf('ang=%f\ti=%f\tj=%f\n',ang,i,j);
            i0 = floor(i);
            j0 = floor(j);
            i1 = i0+1;
            j1 = j0+1;
            %rw = (w2-r)/w2;
            %rw = 1;
            rw = h(round(5*(w2+r)));
            wi0 = (i1-i);
            wi1 = (i-i0);
            wj0 = (j1-j);
            wj1 = (j-j0);
            % wrap along v. axis
            j0 = mod(j0,w);
            j1 = mod(j1,w);
            mask(i0+1,j0+1) = mask(i0+1,j0+1) + wi0*wj0*rw;
            mask(i1+1,j0+1) = mask(i1+1,j0+1) + wi1*wj0*rw;
            mask(i1+1,j1+1) = mask(i1+1,j1+1) + wi1*wj1*rw;
            mask(i0+1,j1+1) = mask(i0+1,j1+1) + wi0*wj1*rw;
        end % traverse a line
        S=sum(mask(:));
        mask = mask*(1/S);
        fourier_mask(:,ai) = mask(:);
    end
end