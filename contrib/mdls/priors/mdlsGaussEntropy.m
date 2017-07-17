%
% Returns the relative entropy (in BITS) of a continuous r.v. with Gaussian distribution
% given its variance.
%
function h=mdlsGaussEntropy(sigma2)
    h=0.5*log2(2*pi*exp(1)*sigma2);
end