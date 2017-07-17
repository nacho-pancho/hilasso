function [Xd,Ao] = compute_ols(Z,D,A,reg,thres)
  N = size(A,2);
  [M,K] = size(D);
  Xd = zeros(M,N);
  Ao = sparse(zeros(size(A)));
  if ~exist('thres','var')
      thres = 0.005;
  end
  e = sum(abs(A))*thres
  if ~exist('reg','var')
      reg = 0;
  end
  for j = 1:N
    idx = find(abs(A(:,j)) > e(j));
    Dg = D(:,idx);
    k = size(Dg,2);
    if length(idx)>0 && length(idx)<M/2
        if reg > 0
            Ao(idx,j) = inv(Dg'*Dg+reg*eye(k)) * Dg' * Z(:,j);
        else
            Ao(idx,j) = inv(Dg'*Dg) * Dg' * Z(:,j);
        end
        Xd(:,j) = Dg * Ao(idx,j);
    else
        Xd(:,j) = D*A(:,j);
    end
  end
end