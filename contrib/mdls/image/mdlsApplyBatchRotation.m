function Y = mdlsApplyBatchRotation(X,RR,angmap)
    J=size(RR,2);
    if ~exist('angmap','var')
        % 
        % mode 1: create J rotated images for each input image
        %
        X = X(:);
        W = sqrt(length(X));
        Y = zeros(W,W,J);
        for j=1:J
            Yj = reshape(RR(:,j),W^2,W^2)*X;
            Y(:,:,j) = reshape(Yj,W,W);
        end
    else
        %
        % mode 2: use angmap to rotate each patch in X
        % and save the rotated in Y
        %
        W = sqrt(size(X,1));
        N = size(X,2);
        Y = zeros(size(X));
        for i=1:N
            j=angmap(i);
            Y(:,i) = reshape(RR(:,j),W^2,W^2)*X(:,i);
        end
    end
end
