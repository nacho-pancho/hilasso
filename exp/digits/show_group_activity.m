function sga=show_group_activity(ga)
  sga = zeros(1,length(ga));
  M =  max([max(ga) 0.1]);
  T1 = 0.8*M;
  T2 = 0.1*M;
  for i=1:length(ga)
      if ga(i) >= T1
          sga(i) = '^';
      elseif ga(i) >= T2 
          sga(i) = '-';
      elseif ga(i) > eps
          sga(i) = '.';
      else
          sga(i) = '_';
  end
  sga = char(sga);
end