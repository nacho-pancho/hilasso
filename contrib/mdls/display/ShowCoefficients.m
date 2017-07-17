function ShowCoefficients(A)
subplot(2,1,1)
imagesc(A)
colormap hot
colorbar
subplot(2,1,2)
[c,h] = hist(A(:),1000);
%bar(h,log(c),'b')
semilogy(h,c); grid on;
