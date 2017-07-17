function Y = rand_sampling(training, num_smp)
% sample local features for unsupervised codebook training

clabel = unique(training.label);
nclass = length(clabel);

for t=1:nclass
    
    num_img = length(find(training.label==clabel(t))); % num of images
    num_per_img = round(num_smp/num_img);
    num_smp = num_per_img*num_img;
    
    idx = find(training.label==clabel(t));
    

    load(training.path{idx(1)});
    dimFea = size(feaSet.feaArr, 1);

    % X = zeros(dimFea, num_smp);
    X = [];
    cnt = 0;

    for ii = 1:num_img,
        fpath = training.path{idx(ii)};
        load(fpath);
        num_fea = size(feaSet.feaArr, 2);
        rndidx = randperm(num_fea);
        num_per_img2 = min(num_per_img,length(rndidx));
        %     X(:, cnt+1:cnt+num_per_img2) = feaSet.feaArr(:,
        %     rndidx(1:num_per_img2));
        X = [X feaSet.feaArr(:, rndidx(1:num_per_img2))];
        cnt = cnt+num_per_img;
    end;
    
    Y{t} = X;
    
end
