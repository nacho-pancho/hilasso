for i=1:size(errorL,2)
    for j=1:size(errorL,3)
        ML(i,j)=min(errorL(:,i,j));
        mL(i,j)=min(mseL(:,i,j));

        MH(i,j)=min(errorH(:,i,j));
        mH(i,j)=min(mseH(:,i,j));
 
        MG(i,j)=min(errorG(:,i,j));
        mG(i,j)=min(mseG(:,i,j));

        aL(i,j)=min(aSetL(:,i,j));
        aH(i,j)=min(aSetH(:,i,j));
        aG(i,j)=min(aSetG(:,i,j));
        
        aux = errorC(:,:,i,j);
        MC(i,j) = min(aux(:));
        
        
        aux = mseC(:,:,i,j);
        mC(i,j) = min(aux(:));
        
        
        aux = aSetC(:,:,i,j);
        aC(i,j) = min(aux(:));
        
        
        
    end
end