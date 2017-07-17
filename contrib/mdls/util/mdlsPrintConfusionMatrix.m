function mdlsPrintConfusionMatrix(cm,labels)
    % header
    NC=length(labels);
    fprintf('tc\\ec\t&');
    for c=1:NC
        fprintf('%4d\t&',labels(c));        
    end
    % body
    fprintf(' e\t\\\\\n');
    for tc=1:NC
        fprintf('%4d\t&',labels(tc));        
        for ec=1:NC
            fprintf('%6.1f &',cm(tc,ec));        
        end        
        fprintf('%5.2f\t\\\\\n',cm(tc,NC+1));
    end
    %
    % footer
    %
    fprintf(' e\t&');
    for ec=1:(NC+1)
        fprintf('%5.2f\t&',cm(tc+1,ec));        
    end
    fprintf('\\\\\n');
end