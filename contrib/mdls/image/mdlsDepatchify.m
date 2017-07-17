%
% function I=depatchify(M,N,X,overlap,DC)
% 
% purpose: 
% reconstruct an image from a set of patches
%
% input:
% M ........ width of image
% N ........ height of image
% X ........ set of square w*w patches (with DC removed)
% overlap ..length of overlapping region
% DC ....... vector of DC coefficients of the patches in X
%
% output:
% I ........ depatchified image
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function I=mdlsDepatchify(M,N,X,overlap,DC)
I=zeros(M,N);
[d,n]=size(X);
add_dc=0;
if nargin > 4
    add_dc=1;
end

width=sqrt(d); % assume square patches!
dw=width-overlap;
%
% overlapping weighting function
%
ow=ones(1,width);
step=1/(overlap+1);
ow(1:overlap)=step*(1:overlap);
ow(width-overlap+1:width)=ow(overlap:-1:1);
wf=ow'* ow;

%n=floor((M-overlap)/dw)*floor((N-overlap)/dw);
% 
% create X from overlappedpatches
%
ii=1;
for i=1:dw:M-width+1
  for j=1:dw:N-width+1
     patch=reshape(X(:,ii),width,width);
    if  add_dc
        patch=patch+DC(ii);
    end
    I(i:i+width-1,j:j+width-1)=I(i:i+width-1,j:j+width-1)+ wf .* patch;        
    ii=ii+1;
  end
end
