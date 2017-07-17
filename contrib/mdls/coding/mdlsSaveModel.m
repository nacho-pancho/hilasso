%
% function modelname=mdlsSaveModel(dirname,imgname,width,overlap,D,metric,lambda)
% 
function modelname=mdlsSaveModel(dirname,imgname,width,overlap,D,alpha,metric,lambda)
  K=size(D,2);
  suffix=mdlsModelSuffix(imgname,width,overlap,K,metric,lambda);
  modelname=sprintf('model_%s.mat',suffix);
  save([dirname '/' modelname],'D','alpha','width','overlap','imgname','lambda','metric');
