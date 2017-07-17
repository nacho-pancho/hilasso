%
% function cm = mdlsConfusionMatrix(label,true_label)
%
% Utility function to compute the Confusion Matrix of a classification
% system.
% 
% Input:
%
% label ............... labels computed by the system
% true_label .......... true labeling
%
% Output:
% 
% cm .................. the confusion matrix
%
function cm = mdlsConfusionMatrix(est_label,true_label)
    class_labels = sort(unique(true_label));
    Nt = length(true_label);
    C = length(class_labels);
    cm = zeros(C,C);
    for t=1:Nt
        tli = find(class_labels == true_label(t));
        eli = find(class_labels == est_label(t));
        cm(tli,eli) = cm(tli,eli)+1;
    end
    conf_row = (100*diag(cm)./(sum(cm,2)+eps));
    conf_col = (100*(diag(cm)')./(sum(cm)+eps));
    conf_tot = (100*trace(cm)/(sum(cm(:))+eps));
    cm = [cm       conf_row; ...
          conf_col conf_tot ];
end