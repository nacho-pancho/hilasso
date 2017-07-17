% function stats = mdlsImageCoeffStatsLearn(D,images_folder[,list_of_images[,stats]])
%
% Specialized version of mdlsCoeffStatsLearn for images
% it will take a dictionary, a list of images and learn 
% statistics of recovery coefficients from all the patches
% in the presented images.
%
% INPUT
%
% D ................. dictionary
% images_folder ..... folder where images are stored. All images are
%                     processed, assuming PGM first and then PNG if no
%                     PGM files are present.
% list_of_images .... optionally, give the explicit list of
%                     images. Useful for debugging.
% stats ............. optional, statistics from previous run
%
% OUTPUT
%
% stats ............. updated statistics
%
%
%======================================================================
% PARSE ARGUMENTS
%======================================================================
%
function stats = mdlsImageCoeffStatsLearn(D,images_folder,list_of_images,stats)
if ~exist('D','var')
    load results/dictionaries/pascal-M0064-k0256-mu0-xmu0-dict.mat
    D = D{1};
    warning('Using default 64x64 256 atom dictionary');
end
if ~exist('images_folder','var')
    images_folder = 'data/images';
    list_of_images = {'boat.png'};
    warning('Running in test mode using boat.png');
end
if ~exist('list_of_images','var')
    list_of_images = mdlsGetFileList(images_folder,'pgm');
    if isempty(list_of_images)
        list_of_images = mdlsGetFileList(images_folder,'png');
        if isempty(list_of_images)
            warning('No images given.');
            return;
        end
    end
end
if ~exist('stats','var')
    stats = mdlsCoeffStatsLearn(D);
end
%
%======================================================================
% INIT
%======================================================================
%

w = 8;
K = size(D,2);
M = size(D,1);
ov = w-1;
%
% maximum rounding error corresponding to 8bpp images
% is about 1/2^(8+1)
% times number of pixels
qe = w^2*(2^-9)^2;
cachedir = 'cache/image_coeffs';
if ~exist(cachedir,'file')
    fprintf('Creating folder %s', cachedir);
    mkdir(cachedir);
end

stats = struct();
stats.AAt = zeros(K,K);
stats.XAt = zeros(M,K);
stats.muA = zeros(K,1);
stats.muA1 = zeros(K,1);
stats.muA2 = zeros(K,1);
stats.AAt = zeros(K,K);
stats.rho = zeros(K,1);
stats.xrho= zeros(K,K);
stats.N   = 0;

%
%======================================================================
% MAIN LOOP
%======================================================================
%
massage_par = mdlsMassageData()
massage_par.normalize = false;
NI=length(list_of_images)
for i=1:NI
    img = list_of_images{i};
    fprintf('%05d/%05d %s ',i,NI,img);
    I = imread(sprintf('%s/%s',images_folder,img));
    
    cachefile = sprintf('%s-%s-w%02d-k%04d-o%02d.mat', ...
                        images_folder,img,w,K,ov);
    cachefile(find(cachefile=='/')) = '-';
    cachefile = [cachedir '/' cachefile];
    X = mdlsDeconstructFast(I,w,ov); % 50% overlap for now. Needs to
                                     % be 100% for learning spatial relations
    X = mdlsMassageData(double(X),massage_par);
    fprintf('N(i)=%d                     \r',size(X,2));
    if ~exist(cachefile,'file')
        N = size(X,2);
        %==================================================================
        % L1 REGULARIZED CODING
        %==================================================================
        %
        s = min(K,M-1);
        l1_mode = 1; % error-constrained
        A = scLasso(X,D,s,qe,l1_mode);
        save(cachefile,'A','N');
    else
        load(cachefile);
    end
    %
    %======================================================================
    % UPDATE STATISTICS
    %======================================================================
    %
if 0 % control to see how much quantization we can afford
    mdlsFigure('Eq'); hold off;
    paletta='ygcmrbk';
    leg={};
    for q=7:-1:3
        Aq=quant(A,2^-q); Eq=X-D*Aq; [h,c] = hist(Eq(:),1000); 
        semilogy(c,h,paletta(q)); hold on;
        leg = {leg{:} sprintf('q=%d',q)};
    end
    legend(leg);
    keyboard
end
    A = quant(A,2^-4); % very little extra distortion
    stats.AAt = stats.AAt + A*A';
    stats.XAt = stats.XAt + X*A';
    stats.muA  = stats.muA + sum(A,2);
    stats.muA1  = stats.muA1 + sum(abs(A),2);    
    stats.muA2  = stats.muA2 + sum(A.^2,2);    
    A = abs(sign(A));
    stats.rho = stats.rho + sum(A,2);
    stats.xrho = stats.xrho + A*A';
    stats.N = stats.N + N;
end
%
%======================================================================
% SHOW THE STUFF
%======================================================================
%
%
% sort atoms by frequency
%
[rhos,idxs] = sort(stats.rho,'descend');
Ds = D(:,idxs);
mdlsFigure('Dictionary');
imagesc(mdlsDictDisplay(Ds));
colormap gray;

mdlsFigure('Frequency (sorted)');
plot(rhos/stats.N);
%
% joint occurence is also rearranged so that most 
% frequent atoms show up first
%
xrhos = stats.xrho(idxs,idxs);
mdlsFigure('Joint frequency');
imagesc(xrhos/stats.N);
colormap cool;

xrhos = stats.xrho(idxs,idxs);
mdlsFigure('Conditional frequency');
xrhos2= xrhos./repmat(diag(xrhos),1,K);
imagesc(xrhos2);
colormap hot;

mdlsFigure('AAt');
imagesc(stats.AAt(idxs,idxs)/stats.N);
colormap hot;
