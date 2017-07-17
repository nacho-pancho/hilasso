function [beta] = ScSPM(feaSet, B, gamma, beta, pyramid, weights, pooling)
%================================================
% 
% Usage:
% Compute the linear spatial pyramid feature using sparse coding. 
%
% Inputss:
% feaSet        -structure defining the feature set of an image   
%                   .feaArr     local feature array extracted from the
%                               image, column-wise
%                   .x          x locations of each local feature, 2nd
%                               dimension of the matrix
%                   .y          y locations of each local feature, 1st
%                               dimension of the matrix
%                   .width      width of the image
%                   .height     height of the image
% B             -sparse dictionary, column-wise
% gamma         -sparsity regularization parameter
% beta          -smooth regularization for sparse coding
% pyramid       -defines structure of pyramid 
% weights       -defines weight for feature on each level
% pooling       -pooling methods for computing the feature
% 
% Output:
% beta          -multiscale sparse coding feature
%
% Written by Jianchao Yang @ NEC Research Lab America (Cupertino)
% Mentor: Kai Yu
% July 2008
%
% Revised Sep. 2009
%===============================================

dSize = size(B, 2);
nSmp = size(feaSet.feaArr, 2);
img_width = feaSet.width;
img_height = feaSet.height;
idxBin = zeros(nSmp, 1);

sparseCodes = zeros(dSize, nSmp);

% compute the sparse code for each local feature
A = B'*B + 2*beta*eye(dSize);
Q = -B'*feaSet.feaArr;

for iter1 = 1:nSmp,
    sparseCodes(:, iter1) = L1QP_FeatureSign_yang(gamma, A, Q(:, iter1));
end;

% compute the pyramid feature
pLevels = length(pyramid);
% bins on each level
pBins = pyramid.^2;
% total bins
tBins = sum(pBins);

beta = zeros(tBins*dSize, 1);

for iter1 = 1:pLevels,
    nBins = pBins(iter1);
    
    pbeta = zeros(nBins*dSize, 1);
    
    wUnit = img_width / pyramid(iter1);
    hUnit = img_height / pyramid(iter1);
    
    % find to which block (bin) each feature belongs
    xBin = ceil(feaSet.x / wUnit);
    yBin = ceil(feaSet.y / hUnit);
    idxBin = (yBin - 1)*pyramid(iter1) + xBin;
    
    for iter2 = 1:nBins,
        sidxBin = find(idxBin == iter2);
        if isempty(sidxBin),
            continue;
        end;
        
        switch pooling,
            case 'energy'    
                % square root of mean of energy pooling
                pbeta((iter2-1)*dSize + (1:dSize)) = sqrt(mean(sparseCodes(:, sidxBin).^2, 2));
            case 'max'
                % max pooling of absolute value
                pbeta((iter2-1)*dSize + (1:dSize)) = max(abs(sparseCodes(:, sidxBin)), [], 2);
            case 'absolute'
                % mean pooing of absolution value
                pbeta((iter2-1)*dSize + (1:dSize)) = mean(abs(sparseCodes(:, sidxBin)), 2);
            otherwise
                disp('unkown pooling method!!!');
                error;
        end;
    end
    
    pbeta = pbeta./sqrt(sum(pbeta.^2));
    beta(sum(pBins(1:iter1-1))*dSize + (1:nBins*dSize)) = pbeta*weights(iter1);   
end;
