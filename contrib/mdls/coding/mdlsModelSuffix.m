function s=mdlsModelSuffix(imgname,width,overlap,K,metric,lambda)
  s=sprintf('%s_L%1d_w%02d_o%02d_k%04d_l%0.5f',imgname,metric,width,overlap,K,lambda);
end