%
% takes various measures on usage
% the usage matrix is received as a set of C columns, where each column has the
% usage of each of the class K by class.
function stats=AnalyzeUsage(usage)
  [K,C]=size(usage);
  mu = sum(usage,2)/C;
  u2=usage-repmat(mu,1,C);
  % covariance matrix: rank C
  Sigma=u2*u2'*(1/(C-1));
  lambda=eig(Sigma);
  lambda=lambda(end-C+1:end);
  stats = struct('eigenvalues',lambda);
  %
  % some geometric measures of scattering
  %
  stats.trace = sum(lambda)/C;
  stats.area = prod(lambda)^(1/C);
  %
  % minimum distance between two vectors
  %
  md = realmax;
  ad = 0;
  ds = [];
  k = 0;
  for i=1:C
    for j=i+1:C
      d=usage(:,i)-usage(:,j);
      d=sqrt(d'*d);
      %d=sum(abs(d)); %sqrt(d'*d);
      ad = ad + d;
      ds = [ds d];
      if d < md
        md = d;
      end
      k = k + 1;
    end
  end
  stats.minDist = md;
  stats.avgDist = ad/k;
  stats.medianDist = median(ds); 
end
