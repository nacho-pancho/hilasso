%
% function [D,A]=dicLearnL1(X,D0,steps,lambda,outdir)
%
% purpose: 
% trains a dictionary and obtains representation coefficients
% for the following sparse problem 
%
% min ||X-D*A||_1 + lambda*||A||_1
% s.t. 
% ||D_i||_2 = 1, for all 1 <= r <= k
%
% The optimization is carried out using alternate minimization
% of the objective function w.r.t. D and A respectively.
% This version implements each optimization step using a fixed point
% technique, but gradient descent or Newton-Rhapson could also be
% used and implemented in future versions.
%
% inputs:
% X ....... data, d x n
% D0 ...... previous dictionary, d x k
% steps ... number of training steps
% lambda .. penalty term for regularization 
% eta ..... quasi-norm smoothing coefficient (small)
% method .. variants on optimization method.
%           This is a bit field with the following structure (defaults to 0x01):
%           bit 1: 1 means eta is updated adaptively. In any case, it always
%                  begins with the value given as argument.
%           bit 2: the same for lambda.
%           bit 3: 1 means use Lasso for sparse coding stage. This is faster and
%                  actually gives sparse results. 0 is consistent with the
%                  L1 prior on error, but gives dense A vectors for some
%                  (yet not understood) reason.
%           bit 4: do not recompute A from scratch at each iteration. Useful
%                  only with my sparse coding stage.
%
% dicdir .. if specified, save dictionary at each step as an
%           image under specified directory.
%
% outputs:
% D ....... final dictionary
% A ... final representation coefficients.
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function [D,A]=mdlsDictLearnL1(X,D0,steps,lambda,eta,method,dicdir)

if ~exist('eta','var')
  eta=1e-5;
end

if ~exist('method','var')
  method=uint8(1);
end


D=D0;
%sparsity=size(X,1);
sparsity=floor(sqrt(size(X,1)));
%sparsity=floor(size(X,1)/2);
% figure(1); colormap(gray); mdlsDicDisplay(D); refresh(1);
%myprint(['.'],sprintf('d%05d',0));
%imwrite(dicDisplay(D),sprintf('d%05d.png',j),'png');

do_diagnostics=exist('dicdir','var');

if do_diagnostics
  if ~exist(dicdir,'file')
    mkdir(dicdir);
  end
  imwrite(mdlsDicDisplay(D,2,255,4),sprintf('%s/d%05d.png',dicdir,i),'png');
end
[d,n]=size(X);
[d,k]=size(D);
%A=zeros(k,n);
A=sign(2*rand(k,n)-1); % stupid initialization!
dd=Inf;
t0=cputime();
error_energy=zeros(1,steps);
A_energy=zeros(1,steps);
total_energy=zeros(1,steps);
laplacianity=zeros(1,steps);
dx=.01;
hcenters=[-1:dx:1];
dD_emp=zeros(length(hcenters),steps);
dD_lap=zeros(length(hcenters),steps);
KL_div=zeros(1,steps);
legs=cell(1,steps);
S=ceil(sqrt(d));
for ss=1:steps
  % 
  % energy from previous step
  %
  E=X-D*A;
  error_energy(ss) = sum(sum(abs(E)));
  A_energy(ss) = sum(sum(abs(A)));
  total_energy(ss) = error_energy(ss) + lambda * A_energy(ss);
  fprintf('%02d: E = %f + %1.4f * %f = %f\n------\n', ...
    ss-1, error_energy(ss), lambda, A_energy(ss), total_energy(ss));
  %
  % 1. sparse coding (A update)
  %
  fprintf('%02d: SPARSE CODING ....... ',ss);
  if bitand(method,4)
    % Lasso: L1 constrain, L2 fitting
    A=mexLasso(X,D,S,lambda,2);
  else
    if bitand(method,8)
      % mine : do not reset A
      t1=cputime();
      A1=A;
      for kk=1:5
        A1=mdlsSparseCodingL1Step(X,D,A1,sparsity,lambda,1e-4,eta);  %thres=0
      end
      E=X-D*A1;God blesse
      A=A1;
    else
      % mine: recompute A
      A=sparseCodingL1(X,D,sparsity,lambda);
    end
    fprintf('dt=%6f |d_A|_1=%5f |E|_1=%5f\n',...
      cputime()-t1,...
      sum(sum(abs(A1-A)))/(k*n), ...
      sum(sum(abs(E)))/(d*n));
  end
  %
  % 2. dictionary update
  %
  fprintf('%02d: DICTIONARY UPDATE ... ',ss);
  t1=cputime();
  for kk=1:5
    D1=mdlsDicLearnL1Step(X,D,A,eta);
  end
  E=X-D1*A;
  fprintf( 'dt=%6f     |d_D|_1=%5f |E|_1=%5f\n', ...
    cputime()-t1,...
    sum(sum(abs(D1-D)))/(d*k), ...
    sum(sum(abs(E)))/(d*n));
  D=D1;
  if bitand(method,1)
    eta=error_energy(ss)/(d*n*100);
  end
  if bitand(method,2)
    lambda=(k*n)/A_energy(ss);
  end
  fprintf('Eta=%f lambda=%f\n',eta,lambda);

  %
  % diagnostic output
  %
  if do_diagnostics
    imwrite(mdlsDicDisplay(D,2,255,4),sprintf('%s/d%05d.png',dicdir,ss),'png');
    %figByName('dic'); colormap(gray); dicDisplay(D); refresh(1);
    dD=unroll(mdlsDicPredict(D));
    
    dD_emp(:,ss)=(hist(dD,hcenters)' +1)/(k*d + 1);
    lambda_ML=mdlsLaplacianFit(dD);
    dD_lap(:,ss)=dx*mdlsLaplacianEval(hcenters,lambda_ML)';
    KL_div(1,ss)=mdlsKLDiv(dD_emp(:,ss)',dD_lap(:,ss)');
    legs{1,ss}=sprintf('%02d:KL=%0.4f',ss,KL_div(1,ss));
  end
  %figure(2); colormap(hot); dicDisplay(DD); refresh(2);
  %
  % 3. stopping condition.
  %
  %myprint(['.'],sprintf('d%05d',j));
end
if do_diagnostics
  figByName('energy');
  plot(1:steps,[total_energy; error_energy; A_energy]);
  legend('Total','Error','Alpha');
  xlabel('iteration');
  title('Dictionary training: evolution of energy term.');
  ylabel('energy');
  fprintf('Total time: %6f\n',cputime()-t0);

  figByName('Laplacianity');
  plot(hcenters,dD_emp); hold on;
  plot(hcenters,dD_lap); hold off;
  title('Laplacianity test');
  legend(legs);
  xlabel('x');
  ylabel('p(x)');
end

