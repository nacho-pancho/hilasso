%%
%% Plots the Receiver Operating Characteristic (ROC) of a 
%% two class hyphothesis testing detector.
%%
%% inputs:
%%
%% The data is given as fp=ROC(fp)
%% false_positives (fp) ... the x axis values 
%% true_positives (tp) .... the y axix values
%%
%% outputs:
%%
%% f ...................... handle to the new figure
%%
function [f]=plotROC(false_positives,true_positives,fprev)
  if (nargin <3)
    f=figure(); 
    xlabel('False Positives');
    ylabel('True Positives');
    title('ROC');
    %hold on, plot([0,1],[0,1],'g'), hold off
  else
    hold on;
    f=fprev;
  end
  plot(false_positives,true_positives,'-b');
  hold off
end