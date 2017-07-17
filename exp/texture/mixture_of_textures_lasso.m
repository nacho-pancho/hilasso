%
% setup
%
N = 100;
M = 100;
img_size = 'smaller';
lambda1v = [ 0.0625];
%lambda2v = [ 0.05 0.1 0.2 ];

NACTIVE  = 2;
NGROUPS  = 8;
NREPS    = 10;
texdir   = ['data/brodatz/' img_size];
textures = mdlsGetFileList(texdir,'pgm','D');

NTEX     = length(textures);
patch_width = 10;
overlap     = patch_width-1;
%aux         = randperm(NTEX);
%chosen      = aux(1:NGROUPS);
%all_sources = textures(chosen);
all_sources = {'D49.pgm',... % 1
               'D84.pgm',... % 2
               'D53.pgm',... % 3
               'D52.pgm',... % 4
               'D33.pgm',... % 5
               'D3.pgm', ... % 6
               'D24.pgm',... % 7
               'D6.pgm'};    % 8
%
% Hand-picked
%
active_idx_v = [ 6 8; ...
                 7 8; ...
                 8 5; ...
                 7 4; ...
                 4 5];
%
% all 14 possibilities
%
active_idx_v = [ 1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8;...
                      2 3; 2 4; 2 5; 2 6; 2 7; 2 8;...
                           3 4; 3 5; 3 6; 3 7; 3 8;...
                                4 5; 4 6; 4 7; 4 8;...
                                     5 6; 5 7; 5 8;...
                                          6 7; 6 8;...
                                               7 8];

active_idx_v = [ 1 1; 1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8;...
                      2 2; 2 3; 2 4; 2 5; 2 6; 2 7; 2 8;...
                           3 3; 3 4; 3 5; 3 6; 3 7; 3 8;...
                                4 4; 4 5; 4 6; 4 7; 4 8;...
                                     5 5; 5 6; 5 7; 5 8;...
                                          6 6; 6 7; 6 8;...
                                               7 7; 7 8];
active_idx_v = [ 1; 2; 3; 4; 5; 6; 7; 8];

NREPS = size(active_idx_v,1);

massage_params = mdlsMassageData();
massage_params.remove_dc = 0;
massage_params.normalize = 0;

rand('twister',8865929);
randn('state' ,8865929);

outdir = sprintf('results/lasso/brodatz/%s-w%02d-ov%02d',img_size,patch_width,overlap);
system(['mkdir -p ' outdir]);

for i=1:NGROUPS
    system(['cp ' texdir '/' all_sources{i} ' ' outdir '/' all_sources{i}]);
end

diary lasso_textures.txt;
%
% load dictionaries for textures that are potentially present
% 
Do    = cell(1,NGROUPS);
mp    = mdlsGetModelPrefix();
DCatom = (1/patch_width)*ones(patch_width^2,1);

for i=1:NGROUPS
    mp.M          = patch_width^2;
    mp.K          = 300;
    texfile       = all_sources{i};
    mp.output_dir = ['results/dictionaries/brodatz/' img_size ];
    mp.base_name  = texfile(1:(end-4));
    prefix        = mdlsGetModelPrefix(mp);
    load([prefix '-dict.mat']);
    Do{i} = D{1};
    %    Do{i} = [ DCatom D{1} ];
end

%
% repetitions
%
fh = fopen('results/mixture_of_textures_lasso.txt','w');
for r = 1:NREPS
    %
    % sample two textures and combine them
    %
    outdir2 = [outdir '/samples'];
    active_idx = unique(active_idx_v(r,:));
    active_sources = all_sources(active_idx);
    NACTIVE = length(active_sources);
    I = cell(1,NACTIVE);
    X = cell(1,NACTIVE);
    DC= cell(1,NACTIVE);
    as0 = zeros(length(all_sources),1); % 'active set' by groups
    as0(active_idx) = 1;
    for i = 1:NACTIVE
        texfile = active_sources{i}
        outdir2 = [outdir2 '-' texfile(1:(end-4))];
        tim     = double( imread([texdir '/' texfile]) )*(1/255);
        I{i}    = tim((end-(M-1)):end,(end-(N-1)):end);
        [X{i},DCt] = mdlsDeconstruct(I{i},patch_width,overlap);
        DC{i} = double(DCt); clear DCt;
        if i == 1
            J = I{i};
        else
            J = J + I{i};
        end
        mdlsFigure('Sources');
        subplot(1,NACTIVE+1,i); imagesc(I{i}); colormap gray; axis off;
    end    
    mkdir(outdir2);
    outdir2 = [outdir2 '/'];

    for i=1:NACTIVE
        texfile = active_sources{i};
        imwrite(mdlsMassageData(I{i},massage_params),[outdir2 '/' ...
                            texfile]);
    end

    mdlsFigure('Sources');
    subplot(1,NACTIVE+1,NACTIVE+1);
    imagesc(J); colormap gray; axis off;
    imwrite(mdlsMassageData(J,massage_params),[outdir2 '/mixture.pgm']);
    %
    % grid
    %
    %
    % deconstruct combined texture
    %
    Y  = mdlsDeconstruct(J,patch_width,overlap);
    c = 1e-3*size(Y,2);

    for lambda1=lambda1v
            outdir3 = sprintf('%s/c-lasso-lambda1%g',...
                              outdir2,lambda1);
            if ~exist(outdir3,'file')
                mkdir(outdir3);
            end
            %
            % collaborative Hilasso
            %
            H = [];
            tol = 1e-3;
            max_iter = 100;
            matfile = sprintf('%s/results.mat', outdir3);
            if ~exist(matfile,'file')
                [Xr,A] = lassoMethod(Y,Do,lambda1);
                save(matfile,'Xr','A');
            else
                load(matfile);
            end
            %            
            % reconstructed images
            %
            for i=1:NGROUPS
                j = find(active_idx == i);
                if isempty(j)
                    Ir{i}   = mdlsReconstruct(Xr{i}, M, N, overlap);
                else
                    Ir{i}   = mdlsReconstruct(Xr{i}, M, N, overlap,DC{j});
                end
                texfile = all_sources{i};
                imwrite( mdlsMassageData(Ir{i},massage_params),...
                         [outdir3 '/' texfile(1:(end-4)) '-rec.png'] );
            end
            %
            % reconstructed mixture
            %
            Yr = Ir{active_idx(1)};
            for i=(1+1):NACTIVE
                Yr = Yr + Ir{active_idx(i)};
            end
            %
            % group active set
            %
            asg = zeros(NGROUPS,size(A,2));
            K = size(Do{1},2);
            for g=1:NGROUPS
                idx = ((g-1)*K+1):g*K;
                asg(g,:) = sum(abs(A(idx,:))) > 0;
            end
            hG = mdlsHammingDistance(asg,repmat(as0,1,size(asg,2)));
            imwrite(mdlsMassageData(Yr,massage_params),[outdir3 '/mixture-rec.png']);
            imwrite(full(A==0),[outdir3 '/spy.png']);
            error = separationError(X,Xr);
            for i=1:NACTIVE
                fprintf(fh,'x%d=%s\t',i,active_sources{i});
            end
            fprintf(fh,'lambda1=%g\tlambda2=%g\tAMSE=%f\thamm=%f\n',...
                    lambda1, lambda2, error, hG);
    end

end
fclose(fh);

diary off
