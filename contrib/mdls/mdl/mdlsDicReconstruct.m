function D=mdlsDicReconstruct(DP)
[d,k]=size(DP);
w=sqrt(d);
D=zeros(d,k);
  for k1=1:k
    ap=DP(:,k1);
    a=zeros(1,d);
    p=0;
    for d1=1:w
      if mod(d1,2) % odd, predict l-r
        for d2=1:w
          idx=(d1-1)*w+d2;
          a(idx)=ap(idx)+p;          
          p=a(idx);
        end
      else % predict r-l
        for d2=w:-1:1
          idx=(d1-1)*w+d2;
          a(idx)=ap(idx)+p;          
          p=a(idx);
        end
      end
    end
    D(:,k1)=a;
  end