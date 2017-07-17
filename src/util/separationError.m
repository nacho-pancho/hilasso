
function error = separationError(X,Xr,H)

if ~exist('H','var')
    missing = 0;
else
    missing = 1;
end


N = size(X{1},2);
numSig = length(X);

% Separation Error
error = 0;
if ~missing
    for i=1:numSig
        dif = (X{i}- Xr{i}).^2;
        error = error + sum(dif(:))/numSig/N;
    end
else
    for i=1:numSig
        dif = ((X{i}- Xr{i}).*H).^2;
        difs = sum(dif)./sum(H);
        error = error + sum(difs(:))/numSig/N;
    end
end
