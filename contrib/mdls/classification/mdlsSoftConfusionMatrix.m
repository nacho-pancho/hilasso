%
% function cm = mdlsSoftConfusionMatrix(label,true_label)
%
% Utility function to compute the Confusion Matrix of a soft classification
% system.
% 
% Input:
%
% est_label ............... soft labels computed by the system.
%                       this is of size Cxlength(true_label)
%                       where each column  j has the probability of
%                       sample being of class i for each class i=1..C
% true_label .......... true labeling
%
% Output:
% 
% cm .................. the confusion matrix
%
function cm = mdlsSoftConfusionMatrix(est_label,true_label)
    class_labels = sort(unique(true_label));
    Nt = length(true_label);
    C = length(class_labels);
    cm = zeros(C,C);
    for c=1:C
        cm(c,:) = sum(est_label(:,true_label == class_labels(c)),2)';
    end
    %
    conf_row = (100*diag(cm)./(sum(cm,2)+eps));
    conf_col = (100*(diag(cm)')./(sum(cm)+eps));
    conf_tot = (100*trace(cm)/(sum(cm(:))+eps));
    cm = [cm       conf_row; ...
          conf_col conf_tot ];
end