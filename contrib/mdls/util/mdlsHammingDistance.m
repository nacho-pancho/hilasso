function h = mdlsHammingDistance(A,B)
    As = A ~= 0; clear A;
    Bs = B ~= 0; clear B;
    h = mean(sum(xor(As,Bs)));
end