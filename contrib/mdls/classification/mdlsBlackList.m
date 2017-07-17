function [blacklist,whitelist] = mdlsBlackList(D,ct)
%
% build concatenation of rival dictionaries
%
    C = length(D);
    Dr = cell(1,C);
    for c=1:C
        Dr{c} = [];
        for c2=[(1:(c-1)) ((c+1):C)]
            Dr{c} = [Dr{c} D{c2}];
        end
    end
    %
    % build black list of atoms that are too coherent 
    %
    blacklist = cell(1,C);
    whitelist = cell(1,C); 
    for c=1:C
        aux = Dr{c}'*D{c};
        % aux (1xKc) = max coherence of each atom in D{c} with
        % an atom from any other dictionary.
        aux = max(aux); 
        Kc = size(D{c},2);
        blacklist{c} = find(aux >= ct);
        whitelist{c} = setdiff(1:Kc,blacklist{c});
        fprintf('Class %d has %d blacklisted atoms\n',c,length(blacklist{c}));
    end        
end
