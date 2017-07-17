function mdata=AnalyzeModel(X,D,A)
  [n,N]=size(X);
  E = X-D*A;
  E = E.*E; % squared
  A=full(A(:));
  %
  % distribution parameters
  %
  [theta,kappa,betta]=BMOLFit(A);
  [dummy,lambda]=BerLapFit(A);
  sigma=sqrt( E*(1/(n*N+1)) );
  %
  % empirical distributions
  %
  [empirical,x]=hist(A,1000);
  empirical = empirical*(1/n/N);
  [empiricalError,xerr] = hist(E(:),1000);
  %
  % fitted distributions
  %
  dx=(x(2)-x(1));
  % these functions give the values of the continuous
  % densities f(x). To get discrete PMFs for comparing to 
  % the histograms, we integrate f(x) in an interval of length
  % dx ( the distance between histogram points.)
  % Here the integral is approximated simply as f(x)*dx.
  fittedLaplacian = BerLapEval(x,theta,lambda)*dx;
  fittedMOL = BMOLEval(x,theta,kappa,betta)*dx;
  %
  % Kullback-Leibler Distance from empirical and fitted
  %
  divLap = Divergence(empirical,fittedLaplacian);
  divMOL = Divergence(empirical,fittedMOL);
  mdata=struct('theta',theta,'kappa',kappa,'beta',betta,...
    'lambda',lambda,...
    'sigma',sigma,...
    'evalPoints',x,...
    'empirical',empirical,...
    'fittedLaplacian',fittedLaplacian,...
    'fittedMOL',fittedMOL,...
    'LaplacianKL',divLap,...
    'MOLKL',divMOL);
end
