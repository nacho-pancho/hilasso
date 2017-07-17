function A = mdlsSparsify(A,L)
    [K,N] = size(A);
    [V,I] = sort(abs(A),1,'descend');
    I = I(1:L,:);
    V = V(1:L,:);
    for j=1:N
        V(:,j) = V(:,j).*sign( A(I(:,j),j) );
    end
    J = repmat(1:N,L,1);
    J = J(:);
    I = I(:);
    A = sparse(I,J,V,K,N);
end