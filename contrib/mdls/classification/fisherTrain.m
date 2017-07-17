%%
%% Train a Fisher's Linear Discriminant using data from two
%% classes in X and Y
%%
%% input:
%%
%% X ...... class 1 data, each column corresponds to the features of a single sample
%% Y ...... class 2 data, idem.
%%
%% output:
%%
%% fld .... linear discriminant 
%% detSw .. det(Sw), Sw being the 'within class scatter matrix'
%%
function [fld,detSw] = fisherTrain(X,Y)

[Dim Nx]=size(X);
[Dim Ny]=size(Y);

mux=(sum(X')/Nx)';
muy=(sum(Y')/Ny)';
mu=(sum(X')+sum(Y'))'/(Nx+Ny);

Mx=X-mux*ones(1,Nx);
My=Y-muy*ones(1,Ny);

Sw=Mx*Mx'+My*My';
Sb=Nx*(mux-mu)*(mux-mu)'+Nx*(muy-mu)*(muy-mu)';

%[e1,lambda] = eig(Sb,Sw);
detSw=det(Sw);
if (detSw==0)
  warning("fisherTrain: det(Sw)==0!");
  Sw=Sw+0.1*ones(size(Sw));
end
[e1,lambda] = eig(inv(Sw)*Sb);

[lambda_max,i]=max(diag(lambda));
fld=e1(:,i);
projMuX=mux'*fld;
projMuY=muy'*fld;
if (projMuY < projMuX) 
  printf("invert axis\n");
  fld=-fld;
end
