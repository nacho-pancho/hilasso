
texdir = 'data/brodatz/tiny';
mkdir('results/dictionaries/brodatz');
mkdir('results/dictionaries/brodatz/tiny');

textures = mdlsGetFileList(texdir,'pgm');
NT = length(textures);
patch_width = 10;
overlap     = patch_width-1;
learn_params = mdlsLearnModel();
learn_params.batch_size = 2048;
learn_params.test_size  = 1024;
learn_params.resume     = 1;
learn_params.output_dir = 'results/dictionaries/brodatz/small';

rand('twister',8865929);
randn('state',8865929);

ini_params = mdlsDictIni();
ini_params.K = 300;
ini_params.method = 'patches';
ini_params.layers = 1;
massage_params = mdlsMassageData();
learn_params.debug = 1;
D = cell(1,NT);
for i=1:NT
    texfile = textures{i};    
    fprintf('Learn dictionary for %s\n',texfile);
    im = double(imread([texdir '/' texfile]));
    im = im(:,:,1);
    mdlsFigure('Current texture','nomargin',true);
    imagesc(im); colormap gray;

    [M,N]  = size(im);
    M1 = ceil(M/2);
    N1 = ceil(N*2/3);
    %
    % upper half for learning
    % which is further divided horizontally in 2/3 for learning and 1/3
    % for testing
    %
    im_train = im(1:(M1+overlap), 1:(N1+overlap)); 
    im_xval  = im(1:(M1+overlap), (N1+1-overlap):end);
    X_train  = mdlsDeconstructFast(im_train, patch_width, overlap);
    X_xval   = mdlsDeconstructFast(im_xval , patch_width, overlap);
    X_train  = mdlsMassageData(X_train,massage_params);
    X_xval   = mdlsMassageData(X_xval,massage_params);
    %
    % initialize dictionary
    %
    ini_params.X = X_train;
    learn_params.D0 = mdlsDictIni(ini_params);
    %
    % initialize dictionary
    %
    learn_params.base_name = texfile(1:end-4); % without extension
    learn_params.training_data = X_train;
    learn_params.testing_data = X_xval;   
    %
    %
    %
    Di = mdlsLearnModel(learn_params);
    Di = Di{1};
    di = mdlsDictDisplay(Di);
    D{i} = Di;
    mdlsFigure('Dictionary','nomargin',true);
    imagesc(di); colormap gray; 
    pause(2);
end