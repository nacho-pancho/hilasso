%
% Returns the relative entropy (in BITS) of a continuous r.v. with Laplacian distribution
% given for the parametrization lambda/2*exp(-lambda|x|)
%
function h=mdlsLapEntropy(lambda)
    h=log2(2*exp(1)/lambda);
end