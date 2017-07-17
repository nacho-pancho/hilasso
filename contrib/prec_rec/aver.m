
function [recall,precision]=aver(rec1,prec1,kind)
% Compute the average precision
% kind  - Indices of the queries to be tested

if nargin == 2
    kind = 1:length(rec1);
end

recall=(5:3:99)/100;  % Recall levels for average
for j=1:length(kind)
    kk=kind(j);
    [rec,prec]=modrec(rec1{kk},prec1{kk});      % Modify recall
    pri(:,j)=interp1(rec,prec,recall,'linear');  % Interpolate
end
for i=1:length(recall)
    jj=find(pri(i,:)>0);
    if length(jj) == 0
        precision(i)=nan;
    else
        precision(i)=sum(pri(i,jj))/length(jj);
    end
end

   
