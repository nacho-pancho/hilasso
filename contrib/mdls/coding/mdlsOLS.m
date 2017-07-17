function [Xd,Ao] = mdlsOLS(Z,D,A,reg)
  N = size(A,2);
  [M,K] = size(D);
  Xd = zeros(M,N);
  if ~exist('thres','var')
      thres = mean(abs(A))*1e-3;% sqrt(eps);
  end
  Ao = sparse(zeros(size(A))); 
  if ~exist('reg','var')
      reg = 0;
  end
  for j = 1:N
    idx = find(abs(A(:,j)) > thres(j));
    Dg = D(:,idx);
    k = size(Dg,2);
    if length(idx)>0 && length(idx)<M
        %
        % least squares of overdetermined problem
        %
        if reg > 0
            Ao(idx,j) = inv(Dg'*Dg+reg*eye(k)) * Dg' * Z(:,j);
        else
            Ao(idx,j) = inv(Dg'*Dg) * Dg' * Z(:,j);
        end
        Xd(:,j) = Dg * Ao(idx,j);
    elseif length(idx) >= M
        Xd(:,j) = D*A(:,j);
    else 
        %
        % solution is 0
        %
        Xd(:,j) = 0;
    end
  end
end