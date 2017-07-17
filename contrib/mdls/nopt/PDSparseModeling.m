%
% function D = PDSparseModeling(params,data,D0)
%
% Name: PDSparseModeling
%
% Category: core function
%
% Description: main sparse modeling loop
%
% Input:
% params ... modeling parameters
% data ..... source image data
% D0 ........ initial dictionary (n x K)
% A ........ coefficients (K x N)
%
% Output:
% D ........ updated dictionary
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Idn$
%
function D = PDSparseModeling(params,data,D)
%
% Dimension of the data vectors
%
n = data.patchWidth^2;
%
% Number of atoms in the dictionary
%
K = size(D,2);
%
% Codename
%
codename = params.codename;
%
% Number of data vectors
%
N = size(data.X,2);
if params.debug
  Figure('Original images')
  I = ReconstructImage(data,'ShowImages',true);
end

iter = 1;
fittingCurve = zeros(params.modelingIters,1);
regCurve = zeros(params.modelingIters,1);
cohCurve = zeros(params.modelingIters,1);
cumCohCurve = zeros(params.modelingIters,1);
maxCohCurve = zeros(params.modelingIters,1);
psnrCurve = zeros(params.modelingIters,1);
psnrCurveOMP = zeros(params.modelingIters,1);
costCurve = zeros(params.modelingIters,1);
normCurve = zeros(params.modelingIters,1);
atomChange = zeros(size(D,2));

D00 = D;
TLOOP=cputime();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (iter <= params.modelingIters)
  LogFile(params,'write',sprintf('---------- Iteration %3d ----------\n',iter));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SPARSE CODING FOR MONITORING: OMP
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tic
  D0 = NormalizeDictionary(D);
  A0 = mexOMP(data.X,D0,10); 
  t = toc;
  data2 = data; 
  data2.X = D0*A0;
  psnr = ReconstructionError(data.X,data2.X);
  psnrCurveOMP(iter) = psnr;
  sparseness = full(Sparseness(A0));
  clear A0;
  LogFile(params,'write',sprintf('OMP : time=%3.2f PSNR: %3.2fdB',t,psnr));
  LogFile(params,'write',sprintf(' Sparseness: %3.1f%%%%\n',sparseness));
  %
  % DEBUG: Show reconstructed image
  %
  if params.debug
    Figure('OMP reconstruction')
    ReconstructImage(data2,'ShowImages',true);
    title(sprintf('OMP: Reconstruction and Error (PSNR=%3.2fdB, sparseness %3.2f%%)',psnr,sparseness))
    drawnow
  end
  clear data2;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SPARSE CODING FOR TRAINIG
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tic
  A = PDSparseCoding(data.X,D,1e-2,params);    
  if params.reviveAtoms
    D = ReviveAtoms(D,A);  
  end
  LogFile(params,'write',sprintf('PDSC: time=%3.2f ',toc));
  psnrIni = ReconstructionError(data.X,D*A);
  psnrCurve(iter) = psnrIni;
  sparseness = full(Sparseness(A));
  LogFile(params,'write',sprintf(' PSNR: %3.2fdB',psnrIni));
  LogFile(params,'write',sprintf(' Sparseness: %3.1f%%%%\n',sparseness));
  data2 = data; data2.X = D*A;
  %
  % DEBUG: Show reconstructed image
  %
  if params.debug
    ReconstructImage(data2,'ShowImages',true);
    Figure('PDSC reconstruction')
    ReconstructImage(data2,'ShowImages',true);
    title(sprintf('PDSC: Reconstruction and Error (PSNR=%3.2fdB, sparseness %3.2f%%)',psnrIni,sparseness))
    drawnow
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % COST FUNCTION EVALUATION
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  LogFile(params,'write',sprintf('; FitTerm; RegTerm CohTerm; NormTerm; CostFunc; MutCohe; CumCohe;\n'));
  [cost,fitTerm,regTerm,cohTerm,normTerm]=CostFunction(data.X,D,A,params);
  [mutualCoherence,cumCoherence] = MutualCoherence(D);
  LogFile(params,'write',sprintf('    ; %6.2f; %6.2f; %6.2f; %6.2f; %7.4f; %8.5f; %6.2f;\n',...
    fitTerm,regTerm,cohTerm,normTerm,cost,max(mutualCoherence),max(cumCoherence)));  
  fittingCurve(iter) = fitTerm;
  regCurve(iter) = regTerm;
  cohCurve(iter)= cohTerm;
  maxCohCurve(iter) = max(mutualCoherence);
  cumCohCurve(iter) = max(cumCoherence);
  normCurve(iter) = normTerm; 
  costCurve(iter) = cost;  
  D2save = DisplayDictionary(D,'MarginSize',4,'Magnification',2,'AtomUsage',full(AtomUsage(A)));
  filename = sprintf('%s/Dictionaries/%s-Dictionary-i%02d.png',params.outdir, codename,iter);
  imwrite(D2save,filename)
  filename = sprintf('%s/Dictionaries/%s-Dictionary-i%02d.mat',params.outdir,codename,iter);
  save(filename,'D','A');
  if params.debug
    Figure('Dictionary')
    DisplayDictionary(D,'Magnification',2,'MarginSize',4,'AtomUsage',full(AtomUsage(A)),'AtomChange',atomChange);
    title(sprintf('Iteration %d; %6.2f; %6.2f; %7.4f; %8.5f; %6.2f;\n',...
    iter,fitTerm,cohTerm,cost,max(mutualCoherence),max(cumCoherence)));  
    drawnow
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % DICTIONARY UPDATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  D0 = D;
  tic;
  D = PDDictionaryUpdate(data,D0,A,params);
  t = toc;  
  for k=1:K
    if ( sum(isnan(D(:,k))) > 0 ) 
    	LogFile(params,'write','Atom contains NAN: reinitializing affected atom.\n');
    	D(:,k)= D00(:,k);
    end
  end

  atomChange = D0 - D;
  atomChange = sum(atomChange.*atomChange);
  fprintf('%d:||D(t+1)-D(t)||/(n*K)=%f\n',iter,sum(atomChange)/K/n);
  atomChange = atomChange/max(atomChange);
  %
  % optional forced normalization of dictionary
  %
  if params.forceNormalization > 0
    if ~mod(iter,params.forceNormalization)
      LogFile(params,'write',sprintf('%4d: Normalizing dictionary\n',iter));
      D = NormalizeDictionary(D);
    end
  end
  iter = iter + 1;
end % while, big loop
D2save = DisplayDictionary(D,'MarginSize',4,'Magnification',3,'Plain',1);
filename = sprintf('%s/%s-Final-Dictionary.png',params.outdir,codename);
imwrite(D2save,filename);
if params.debug
  ReconstructImage(data2,'ShowImages',true,'SaveImages',sprintf('%s/RecImg/%s-RecImg-i%02d',params.outdir,codename,iter));
end
hfc = Figure('Curves');
filename = sprintf('%s/%s-data.mat',params.outdir,params.codename);
save(filename,'fittingCurve','regCurve','cohCurve','normCurve','costCurve','psnrCurve','cumCohCurve','maxCohCurve');
