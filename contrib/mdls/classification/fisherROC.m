function [accuracy,FA50,fp,tp] = fisherROC(fld,X,Y)
%
% Receiver Operator Characteristic plot
%
% inputs:
%
% fld ..... Fisher LDA Projection vector
% X ....... Original (cover) features
% Y ....... Stego features
% display . input 0 if no figure with ROC should be displayed, other value otherwise
%
% outputs:
%
% accuracy . area between the ROC curve and the diagonal normalized so that a perfect detection has accuracy=1
% FA50 ..... false alarm (between zero and 1) at 50% of true positives
% fp, tp ... x and y data for the ROC function
%

%
% fld must be a row vector (1xD)
%
if size(fld,1) > size(fld,2) 
  fld=fld';
end

%
% data muse be DxN (any N)
%
if size(X,1) != size(fld,2)
  if size(X,2) == size(fld,2)
    X = X';
  else
    error('fisherROC: X has wrong dimensions');
  end
end

if size(Y,1) != size(fld,2)
  if size(Y,2) == size(fld,2)
    Y = Y';
  else
    error('fisherROC: X has wrong dimensions');
  end
end

%
% X and Y must have 
%
projX =fld*X;
projY=fld*Y;


projX=projX(:)'; projY=projY(:)'; 
leo=length(projX); les=length(projY);
True=[zeros(1,leo),ones(1,les)];   % 0 for original data, 1 for stego data
[alldata,order]=sort([projX,projY]);
N=length(alldata);
TrueSort=True(order);
difference=alldata(2:end)-alldata(1:end-1); 
when_equal=(difference==0); 

j=2;
if mean(projX)<mean(projY)
    ROCx(1)=0; ROCy(1)=0; x=0; y=0;      % false positives and true positives when the threshold is below all values
    when_equal=[0,when_equal];
    for i=1:N
        x=x+(1-TrueSort(N-i+1));
        y=y+TrueSort(N-i+1);
        if ~when_equal(N-i+1), ROCx(j)=x; ROCy(j)=y; j=j+1; end
    end
else
    ROCx(1)=0; ROCy(1)=0; x=0; y=0;      % false positives and true positives when the threshold is below all values
    when_equal=[when_equal,0];
    for i=1:N
        x=x+(1-TrueSort(i));
        y=y+TrueSort(i);
        if ~when_equal(i), ROCx(j)=x; ROCy(j)=y; j=j+1; end
    end
end

fp = ROCx/leo; 
tp = ROCy/les;

accuracy = 2*sum((fp(2:end)-fp(1:end-1)).*(tp(2:end)+tp(1:end-1))/2)-1;

ind_more_than_50 = find(tp>0.5);  ind_more_than_50 = ind_more_than_50(1);
if ind_more_than_50-1<1, 
  Fa50=0; 
end
x1 = fp(ind_more_than_50-1); x2 = fp(ind_more_than_50);
y1 = tp(ind_more_than_50-1); y2 = tp(ind_more_than_50);
FA50 = x1 + (x2-x1)*(0.5-y1)/(y2-y1);

endfunction