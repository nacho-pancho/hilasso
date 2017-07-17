function ge=group_energy(A,gs)
    gs=gs(:);
    if isscalar(gs)
        ge = zeros(size(A,1)/gs,size(A,2));
        for i=1:(size(A,1)/gs)
            idx = (1+(i-1)*gs):(i*gs);
            ge(i,:) = sum( A(idx,:).^2 );
        end
    else        
        ge = zeros(length(unique(gs)),size(A,2));
        for i=1:length(unique(gs))
            idx = find(gs==i);
            ge(i,:) = sum( A(idx,:).^2 );
        end
    end
end
