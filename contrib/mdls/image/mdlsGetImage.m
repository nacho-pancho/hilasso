%
% Get an image from the test database from its canonical
% name (no extension or path) and normalize it
%
function I=mdlsGetImage(canonical_name,normalize)
    if ~exist('normalize','var')
        normalize=true;
    end
    if normalize
        I=double(imread(['data/images/' canonical_name '.png']))* ...
          (1/256);
    else
        I=double(imread(['data/images/' canonical_name '.png']));
    end
end