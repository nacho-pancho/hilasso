%
% function D = mdlsLearnImageDict(params)
%
% Category: image processing
%
% Purpose: General purpose reconstructive dictionary learning for images images.
%
% Description:
%
% Input:
% params ...... algorithm parameters in a single struct.
%               See help mdlsLearModel to see a description of them.
%               This function sets some defaults of its own on these
%               parameters, and adds a few more. These are:
%                 w (new) ...... width of the image patches. Defaults to 8.
%                 K ............ Size of dictionary. Defaults to 256.
%                 batch_size ... Set to 2560 by default.
%                 test_size .... Set to 5000 by default.                
%                 control_step . Set to 5 by default.
%                 normalize .... False. Keep patches as they are.
%                 D0 ........... Empty by default. It is initialized  
%                                as an average of 5 random patches plus
%                                some slight noise.
%                                  corresponding to the data.
% Output:
% D .......... The final dictionary or set of dictionaries. If function
%              is called with no arguments, this is a set of default parameters.
%
function D = mdlsLearnImageDict(params)
    %
    % SET DEFAULTS
    %
    warning(['This function is OBSOLETE and may not work well. Use ' ...
             'mdlsImageDictLearn instead.']);
    %
    % model parameters
    %
    %
    % PARSE INPUT
    %
    if nargin == 0
        % no input means request default parameters
        params = mdlsLearnModel();
        params.w            = 8;
        params.K            = params.w^2*4;
        params.output_dir    = 'results/dictionaries';
        params.training_format    = 'cache';
        params.testing_format    = 'cache';
        params.batch_size = params.K*10;
        params.test_size  = 1024;
        params.D0 = {};
        D = params;
        return;
    end

    if ~isstruct(params)
        error('Input must be a struct. Call function with no arguments for a default set of parameters.');
    end
    %
    % Create training buffer
    %
    training_dir = 'data/pascal06/train';
    training_list = mdlsGetFileList(training_dir,'pgm');
    system('mkdir cache/patches');
    fname = sprintf('cache/patches/training_w%02d-info.mat',params.w);
    if ~exist(fname,'file')
        fprintf('Creating patch buffer (this can take some time)\n');
        compact_mode = 1;
        remove_dc = 0;
        train_buffer_info = mdlsCreatePatchBuffer(training_dir, ...
                                                  training_list,...
                                                  params.w,...
                                                  'cache/patches/training',...
                                                  remove_dc,...
                                                  compact_mode);
    else
        fprintf('Loading train buffer\n');
        load(fname);
        train_buffer_info = buffer_info;
        clear buffer_info;
    end
    %
    % Create testing buffer
    %
    testing_dir = 'data/pascal06/test';
    testing_list = mdlsGetFileList(testing_dir,'pgm');
    fname = sprintf('cache/patches/testing_w%02d-info.mat',params.w);
    if ~exist(fname,'file')
        fprintf('Creating patch buffer (this can take some time)\n');
        compact_mode = 1;
        remove_dc = 0;
        test_buffer_info = mdlsCreatePatchBuffer(testing_dir,...
                                                 testing_list,...
                                                 params.w,...
                                                 'cache/patches/testing',...
                                                 remove_dc,...
                                                 compact_mode);
    else
        fprintf('Loading test buffer\n');
        load(fname);
        test_buffer_info = buffer_info;
        clear buffer_info;
    end
    params.training_data = train_buffer_info;
    params.testing_data  = test_buffer_info;
    %
    % Initial dictionary
    %
    Nb = train_buffer_info.patch_count;
    M  = params.w^2;
    indexes = ceil( rand(1, params.K*2) * Nb );
    indexes=sort(indexes);
    fname = [train_buffer_info.buffer_prefix '-patches.bin'];
    fbuf = fopen(fname);
    if fbuf < 0
        error(sprintf('Error opening %s\n',fname));
    end
    Xb = mdlsGetPatchesFromBuffer(fbuf, train_buffer_info, ...
                                  indexes)*(1.0/255.0);
    Xb = Xb - repmat(mean(Xb),M,1);
    fclose(fbuf);
    
    if isempty(params.D0)
        fprintf('Initializing dictionary\n');
        di_params = mdlsDictIni();
        di_params.K = params.K;
        di_params.dim = size(Xb,1);
        di_params.X = Xb;
        di_params.method = 'patches';
        di_params.sigma = 10/255;
        di_params.layers = 5;
        params.D0 = mdlsDictIni(di_params);
    else
        fprintf('Using a given initial dictionary\n');
    end
    fprintf('Learning dictionary\n');
    D = mdlsLearnModel(params);
end
    
