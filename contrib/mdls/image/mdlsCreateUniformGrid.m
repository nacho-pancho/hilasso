% function grid=mdlsCreateUniformGrid(m,n,dx,dy,x0,y0)
%
% Create sampling grid. For using mainly with mdlsGetPaches.
%
% INPUT
% x0 ...... x coordinate for top left corner point
% y0 ...... y coordinate for top left corner point
% xf ...... x coordinate for bottom right corner point
% yf ...... y coordinate for bottom right corner point
% dx ...... x step
% dy ...... y step
%
function grid=mdlsCreateUniformGrid(x0,y0,xf,yf,dx,dy)
    if ~exist('x0','var')
        x0 = 0;
    end
    if ~exist('y0','var')
        y0 = 0;
    end
    base_x = x0:dx:xf;
    base_y = y0:dy:yf;
    lx = length(base_x);
    ly = length(base_y);
    grid = zeros(2,lx*ly);
    % for a dense grid on a 4x4 image 
    % base_x = [ 1 2 3 4 ] 
    % x = [ 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 ]
    % y = [ 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 ]
    % this is the trick to repeat each element lx times
    % uses fact that data is stored in column major order
    x = repmat(base_x,1,ly);
    grid(1,:) = x; clear x;
    y = repmat(base_y,lx,1);
    grid(2,:) = y(:)';
    grid = int32(grid);
end