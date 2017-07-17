function [re,pr]=modrec(rec,prec)
% Modify the vector REC so that any equal elements are removed
% For equal elements take the average of the PREC elements
j=0; i=0;
while i<length(rec)
    i=i+1;
    j=j+1;
    if i<length(rec) & rec(i+1)<=rec(i) 
        ii=find(rec==rec(i));
        re(j)=rec(i);
        pr(j)=sum(prec(ii))/length(ii);
        i=i+length(ii)-1;
    else
        re(j)=rec(i);
        pr(j)=prec(i);
    end
end
