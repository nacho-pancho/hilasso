function w = simplex_proj(v,z)
  if nargin < 2
      z = 1;
  end
  w = sort(v,'descend');
  t = 1;
  u = w(t)-z;
  %fprintf('t=%d w(t)=%f u=%f\n',t,w(t),u);
  while (t*w(t)) > u
      t = t + 1;
      if t > length(w)
          break;
      end
      u = u + w(t) ;
      %fprintf('t=%d w(t)=%f u=%f\n',t,w(t),u);
  end
  if (t < length(w))
      u = u - w(t);
      t = t - 1;
  end
  theta = u/t;
  w = v - theta;
  w = w.*(w > 0);
end