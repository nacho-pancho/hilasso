

function Dn = getDictionaries(im,lambda,T,numD,sizeD)

sP = 8;
overlap = 7;

[X,Dc] = mdlsPatchify(im,sP,overlap);


Dini = mdlsGenPatchesDictionary(X,256,5/255,10);
Do =learnDict(X,[],{Dini},300,lambda,0,4000);

[S,D] = mdlsInitializeClusters(X,Do{1},numD,T,lambda,sizeD,true);


save intermediate D Do

[clustering Dn Irr] = Texcluster(im,D,T,overlap,0.1,4);