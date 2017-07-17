%
% Get an image from the test database from its canonical
% name (no extension or path) and normalize it
%
function I=GetImage(canonical_name)
    I=double(imread('data/images/' canonical_name '.png'))*(1/255);
end