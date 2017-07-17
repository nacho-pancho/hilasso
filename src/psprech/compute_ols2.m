function [Xd,Ao] = compute_ols2(Z,D,A,reg,thres)
  N = size(A,2);
  [M,K] = size(D);
  Xd = zeros(M,N);
  Ao = sparse(zeros(size(A)));
  if ~exist('thres','var')
      thres = 0.01;
  end
  %nza = A(find(A));
  e = thres; % e = mean(abs(nza))*thres
  if ~exist('reg','var')
      reg = 1e-2;
  end
  for j = 1:N
    idx = find(abs(A(:,j)) > e);
    Dg = D(:,idx);
    k = size(Dg,2);
    if length(idx)>0 % && (length(idx) <= ceil(M/2))        
        if length(idx) > ceil(M/2) % solution is too dense: regularize
                                   % inverse
            fprintf('x');
            Ao(idx,j) = inv(Dg'*Dg+reg*eye(k)) * Dg' * Z(:,j);
        else
            fprintf('.');
            Ao(idx,j) = inv(Dg'*Dg) * Dg' * Z(:,j);
        end
        Xd(:,j) = Dg * Ao(idx,j);
    else 
        fprintf('o');
        Xd(:,j) = D*A(:,j);
    end
  end
  fprintf('\n');
end