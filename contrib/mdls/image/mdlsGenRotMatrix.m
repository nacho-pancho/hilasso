
% parametros
n = 9;
theta = -90;

% inicializo las matrices
B = zeros(n,n,n^2);     % matrices base canonica
Brot = zeros(n,n,n^2);  % matrices base canonica rotadas theta
A = zeros(n^2,n^2);     % matriz A(theta)

% genero las matrices de la base canonica y aprovecho el loop para
% rotarla y agregar la respectiva columna de A(theta)
for i =1:n
    for j =1:n
        
        % creo la matriz
        B(i,j,n*(i-1)+j) = 1;
        
        % roto
        Brot(:,:,n*(i-1)+j) = imrotate(B(:,:,n*(i-1)+j),theta,'bilinear','crop');
        BrotS = reshape(Brot(:,:,n*(i-1)+j),n,n);
        % agrego la columna
        Atheta(:,sub2ind(size(BrotS),i,j))= BrotS(:);
        
    end
end


% testeo que es lo mismo, dada una imagen nxn aplicar imrotate (con
% bilineal) o A(theta)*x

% para cualquier n
% Img = round(rand(n,n));
% para n impar
Img = zeros(n,n); Img(1:(n-1)/2,1:(n-1)/2) = ones((n-1)/2,(n-1)/2);

ImgR1 = imrotate(Img,theta,'bilinear','crop');

ImgR2v = Atheta*Img(:);
ImgR2 = reshape(ImgR2v,n,n);

subplot(1,3,1)
imagesc(Img);colormap gray
axis equal;axis off
title('original')
%figure
subplot(1,3,2)
imagesc(ImgR1);colormap gray
title('imrotate')
axis equal;axis off
subplot(1,3,3)
imagesc(ImgR2);colormap gray
axis equal;axis off
title('A(theta)*x')

sprintf('Diferencia entre las rotaciones: %d',norm(ImgR1-ImgR2))