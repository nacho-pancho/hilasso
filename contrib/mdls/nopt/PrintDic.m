if ~exist('K','var')
  K=64;
end
if ~exist('zeta','var')
  zeta=10;
end

mod=load(sprintf('../output/exp2/MNIST784-%03d-Lasso-sp005-mu0.000-eta00-var0.0100-model.mat',K));
mocod=load(sprintf('../output/exp2/MNIST784-%03d-Lasso-sp005-mu%d.000-eta50-var0.0100-model.mat',K,zeta));

maxD = max(max([mod.D mocod.D]));
minD = min(min([mod.D mocod.D]));
%
% normalize so that both are in the same range between 0 and 1
%
dicMOD   = (mod.D - minD)*(1/(maxD-minD));
dicMOCOD = (mocod.D - minD)*(1/(maxD-minD));
%
% create images
%
iMOD   = DisplayDictionary(dicMOD,'Plain',1);
iMOCOD = DisplayDictionary(dicMOCOD,'Plain',1);
%
% invert them so that background is white
%
iMOD = 255-double(iMOD(:,:,1));
iMOCOD = 255-double(iMOCOD(:,:,1));
%
%
%
imwrite(uint8(iMOD),sprintf('mnist_mod_%03d.png',K));
imwrite(uint8(iMOCOD),sprintf('mnist_mocod_%03d.png',K));
%
% try to equalize them automatically
%
if 0
maxD = max(max([mod.D mocod.D]));
minD = min(min([mod.D mocod.D]));
mod.D   = (mod.D - minD)*(1/(maxD-minD));
mocod.D = (mocod.D - minD)*(1/(maxD-minD));
total = [mod.D mocod.D];

mu=mean(total(:))
sig=var(total(:))
maxD = min(mu+2*sqrt(sig),1)
minD = max(mu-2*sqrt(sig),0)
dicMOD   = (mod.D - minD)*(1/(maxD-minD));
dicMOCOD = (mocod.D - minD)*(1/(maxD-minD));
dicMOD(dicMOD > maxD) = maxD;
dicMOCOD(dicMOCOD > maxD) = maxD;
dicMOD(dicMOD < minD) = minD;
dicMOCOD(dicMOD > minD) = minD;

iMOD   = DisplayDictionary(dicMOD,'Plain',1);
iMOCOD = DisplayDictionary(dicMOCOD,'Plain',1);

iMOD = 255-double(iMOD(:,:,1));
iMOCOD = 255-double(iMOCOD(:,:,1));

imwrite(uint8(iMOD),'mnist_mod_64_auto.png');
imwrite(uint8(iMOCOD),'mnist_mocod_64_auto.png');
end
