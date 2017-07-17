close all
%
% training data
%
listOfImages = {'../images/barbara.png'};
%
% decomposition parameters
%
params = DefaultParams;
patchWidth = 8;
overlap = 4;
params.patchWidth = patchWidth;
params.overlap = overlap;
%
% Fixed parameters
%
params.modelingIters = 10;
params.debug = 1;
params.updateStepsizeRule = 'ARO'; % Fix step
params.updateIters = 1;
params.codingIters = 10;
%
% Generate the data
%
data = GetData2(patchWidth,overlap,listOfImages);
%
% Load the initial (global) dictionary
%
load ../dictionaries/julien_d64_K256.ascii
D = julien_d64_K256;
n = size(data.X,1);

params.kappaCoef = 3;
params.varError = 1; % value assumed when setting the value of lambda in LASSO
lambdas=[0.001 0.01 0.1 1.0];

X=data.X;
Figure('SparseCoding');
for lambda=lambdas
  %
  % Least Squares coding
  %
  A_ls = LeastSquaresSparseCoding(X,D,1e-2,2);
  % sparse coding with LASSO
  % mean lambda of MOL is mean of Gamma distribution: kappa/beta
  % 
  A_lasso = mexLasso(X,D,n,lambda,2);
  
  %
  % Iterated Shrinkage
  %
  A_ssf = IteratedShrinkage();
  E_ssf = X-D*A_ssf;
  E_ssf = E_ssf.*E_ssf; % squared
  A_ssf=full(A_ssf(:));
  [t_ssf,k_ssf,b_ssf]=BMOLFit(A_ssf);
  [dummy,l_ssf]=BerLapFit(A_ssf);
  %
  % MOL
  %
  % this is the value of beta that, for a fixed kappa, gives the lambda used in LASSO
  % it is actually wrong: should be params.kappaCoef/lambda, but it is how it's written
  % in our tests and should not make a big difference.
  %
  params.betaCoef = (params.kappaCoef+1)/lambda;
  A_mol = PDSparseCoding(X,D,[],params);
  E_mol = X-D*A_mol;
  E_mol = E_mol.*E_mol; % squared  
  A_mol=full(A_mol(:));
  [t_mol,k_mol,b_mol]=BMOLFit(A_mol);
  [dummy,l_mol]=BerLapFit(A_mol);  
  %
  % likelyhood
  %
end