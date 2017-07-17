

sigma = [100];%:0.1:0.4;
rep = 3;


lambda1 = [0.1];
lambda2 = [5];
missing = [0.6];

errorC = zeros(length(missing),length(lambda1),length(lambda2));

aSetL   = zeros(length(lambda1));
aSetCGL = zeros(length(lambda2));
aSetCHL = zeros(length(lambda1),length(lambda2));
aSetHL  = zeros(length(lambda1),length(lambda2));

rand('twister',12345888);
%system('mkdir -p results/missing/');
outdir0 = 'results/missing';
load results/missing/dict_pablo.mat
load data/usps/usps_test.mat

active0 = [4 6;...
           1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8; 1 9; 1 10;...
               2 3; 2 4; 2 5; 2 6; 2 7; 2 8; 2 9; 2 10;...
                    3 4; 3 5; 3 6; 3 7; 3 8; 3 9; 3 10;...
                         4 5; 4 6; 4 7; 4 8; 4 9; 4 10;...
                              5 6; 5 7; 5 8; 5 9; 5 10;...
                                   6 7; 6 8; 6 9; 6 10;...
                                        7 8; 7 9; 7 10;...
                                             8 9; 8 10;...
                                                  9 10];
active0 = active0-1;

N = 200;
gs = 300;

for r=1:size(active0,1)
  active = active0(r,:);
  outdir = sprintf('%s/active-%d-%d/',outdir0,active(1),active(2))
  mkdir(outdir);
  as0 = zeros(10,N);
  as0(active+1,:) = 1;

  [Y,X] = createDataDigits(data,N,active);
  H     = rand(size(Y)) < 0.6;
  save([outdir '/data.mat'],'X','Y','H');
  %
  % LASSO
  %
  if 1
      for h=1:length(lambda1)
          fname = sprintf('%s/lasso-lambda1=%g.mat',...
                          outdir,lambda1(h));
          if ~exist(fname,'file')
              [Xr,A,v] = lassoMethod(Y,D,lambda1(h),H);
              save(fname,'Xr','A','v');
          else
              load(fname);
          end
          as=group_act_set(A,gs);
          ge=group_energy(A,gs);
          imwrite(ge,[fname(1:(end-3)) '-group_energy.png']);
          imwrite(double(abs(A)>0),[fname(1:(end-3)) '-act_set.png']);
          aSetL = avg_group_hamming(as,as0);
          fprintf('Lasso: lambda=0.1 hamming=%5.3f seperr=%g\n',aSetL,separationError(X,Xr));
      end
  end

  
  %
  %C-HiLasso
  %
  if 1
      for h = 1:length(lambda1)
          for f=1:length(lambda2)
              fname = sprintf('%s/c_hilasso-lambda1=%g-lambda2=%g.mat',...
                              outdir,lambda1(h),lambda2(f));
              if ~exist(fname,'file')
                  % collaborative lasso method
                  [Xr,A] = HIlassoColMethod(Y,D,lambda1(h),lambda2(f),1e-3,H,200);;
                  save(fname,'Xr','A');
              else
                  load(fname);
              end
              %        errorCHL(m,h,f) = separationError(X,Xr);
              as=group_act_set(A,gs);
              ge=group_energy(A,gs);
              imwrite(ge,[fname(1:(end-3)) '-group_energy.png']);
              imwrite(double(abs(A)>0),[fname(1:(end-3)) '-act_set.png']);
              aSetCHL = avg_group_hamming(as,as0);% mean(sum(AC~=0));
              fprintf('C-HiLasso: lambda1=%g lambda2=%g hamming=%5.3f seperr=%g\n',...
                      lambda1(h),lambda2(f),aSetCHL,separationError(X,Xr));
          end
      end
  end

  %load results/missing/c_hilasso_pablo.mat
  as=group_act_set(A,gs);
  aSetCHL = avg_group_hamming(as,as0);% mean(sum(AC~=0));
                                      %
  fprintf('C-HiLasso: hamming=%5.3f\n',aSetCHL);
  % C-GLasso
  %
  for f=1:length(lambda2)
      fname = sprintf('%s/c_glasso-lambda2=%g.mat',...
                      outdir,lambda2(f));
      if ~exist(fname,'file')
          % collaborative lasso method
          [Xr,A,v]  = GLassoColMethod(Y,D,lambda2(f),1e-6,500);
          save(fname,'Xr','A');
      else
          load(fname);
      end
      %errorGL(m,h,f) = separationError(X,Xr);
      as=group_act_set(A,gs);
      ge=group_energy(A,gs);
      imwrite(ge,[fname(1:(end-3)) '-group_energy.png']);
      imwrite(double(abs(A)>0),[fname(1:(end-3)) '-act_set.png']);
      aSetGL(f) = avg_group_hamming(as,as0); % mean(sum(AC~=0));
      fprintf('C-GLasso: lambda2=%g hamming=%5.3f seperr=%g\n',lambda2(f),aSetGL(f),separationError(X,Xr));
  end

  %
  % HiLasso
  %
  if 0
      for h=1:0% length(lambda1)
          for f=1:length(lambda2)
              fname = sprintf('%s/hilasso-lambda1=%g-lambda2=%g.mat',...
                              outdir,lambda1(h),lambda2(f));
              if ~exist(fname,'file')
                  [Xr,A,v] = HIlassoMethod(Y,D,lambda1(h),lambda2(f));
                  save(fname,'Xr','A');
              else
                  load(fname);
              end
              %errorHL(m,h) = separationError(X,Xr);
              as=group_act_set(A,gs);
              ge=group_energy(A,gs);
              imwrite(ge,[fname(1:(end-3)) '-group_energy.png']);
              imwrite(double(abs(A)>0),[fname(1:(end-3)) '-act_set.png']);
              aSetHL = avg_group_hamming(as,as0); % mean(sum(AC~=0));
              fprintf('C-HiLasso: lambda1=%g lambda2=%g hamming=%5.3f seperr=%g\n',...
                      lambda1(h),lambda2(f),aSetHL,separationError(X,Xr));
          end
      end
  end

end % each set of active numbers



