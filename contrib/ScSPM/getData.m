
direc = 'results/graz/bikes/';
% direc_img = 'data/inria_graz02/bikes/';


nMax = 300;
d = dir([direc '*32.mat']);
% d = d(1:2:nMax);

% datosBg = 100000;

datosBike = 1500000;

Xbike = zeros(128,datosBike);
% Xbg = zeros(128,datosBike);
k=1;
k2=1;
maxImSize = 500;

length(d)


for i = 1:2:nMax
    
    file = [direc d(i).name];
    load(file)

% Codigo que selecciona el Bounding Box
%     mask_file = [direc_img d(i).name(1:end-11) 'matmask.png'];
%     I = imread(mask_file);
%     I = Im{i};
%     
%     [m1,ind1] = max(I);
%     [m2,ind2] = max(I,[],2);
%     
%     id1 = find(m1~=0);
%     id2 = find(m2~=0);
%     
%     I(min(id2):max(id2),min(id1):max(id1)) = 1;
%      
%     for j=1:length(feaSet.x)
%         idx(j) = I(round(feaSet.y(j)),round(feaSet.x(j)))~=0;
%     end
% 
%     ii = find(idx == 1 );
%     L = length(ii);
    
    L = size(feaSet.feaArr,2);
    Xbike(:,k:k+L-1) = feaSet.feaArr;

%     jj = find(idx == 0 );
%     L2 = length(jj);

%     Xbg(:,k2:k2+L2-1) = feaSet.feaArr(:,jj);
%     k2 = k2+L2;

%     aux = randperm(L2);
%     L = max(L2,4000);
%     Xbg(:,k:k+L-1) = feaSet.feaArr(:,jj(aux(1:L)));


    k = k+L;
    disp(L)
    disp(L2)
    disp('---')


    
end

Xbike = Xbike(:,1:k);
% Xbg = Xbg(:,1:k2);
