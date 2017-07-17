%
% Given a patch width W and a set of J different rotation angles, this 
% function returns a W^4 x J sparse matrix where for each J, 
% each column has the weights associated to a single pixel in the original patch
% in the rotated image.
% To obtain a rotated patch Y with angle rotation_angles(j), we do
% Y = RR(:,:,j)*X
%
% INPUT
%
% w ............. patch width
% angles ........ chosen rotation angles, 360
% interp_type ... kind of interpolation done
% 
function RR=mdlsCreateBatchRotation(w,...
                                    rotation_angles, ...
                                    interp_type)
    if ~exist('interp_type','var')
        interp_type = 'bilinear';
    end
    %
    %  see if we don't have it in the cache
    %  angles are assumed to be evenly spaced from first to last
    cachefile = sprintf('cache/rotations/rot-w%02d-%s-angles%g-%g-%g.mat',...
                        w,interp_type,...
                        rotation_angles(1),...
                        rotation_angles(2)-rotation_angles(1),...
                        rotation_angles(end));
    if exist(cachefile,'file')
        load(cachefile);
        return;
    end
    %
    % number of non-zero entries in each matrix
    %
    switch interp_type
      case 'bilinear'
        col_nnz = 4;
      case 'bicubic'
        col_nnz = 16;
    end
    N = w^2;
    M = N^2; % size of each rotation matrix in linear space
    J = length(rotation_angles);
    col_nnz = min(col_nnz,N);
    RR=spalloc(M,J,J*col_nnz);
    %
    % fill the matrix
    %
    A   = zeros(w,w);
    for j=1:J
        markers=['|','/','-','\'];
        for n=1:N
            % A contains a 1 in pixel n (linear indexing)
            A(n)=1;
            % B is the transformation of this pixel
            B = imrotate(A,rotation_angles(j),interp_type,'crop');
            A(n)=0;
            % which we copy in RR, in the n-th column of the j-th trans
            % matrix
            RR( N*(n-1)+find(B), j ) = B(find(B));
        end
        fprintf('\b');
    end
    save(cachefile,'RR','rotation_angles','w');
end