function as=group_act_set(A,gs,thres)
    if ~exist('thres','var')
        thres = sqrt(eps);
    end
    as = zeros(size(A,1)/gs,size(A,2));
    for i=1:(size(A,1)/gs)
        idx = (1+(i-1)*gs):(i*gs);
        as(i,:) = sum( abs(A(idx,:)) ) > thres;
    end
end
