function [prec,tpr, fpr, pred_thresh] = prec_rec(pred,pos,count,num_thresh,varargin)
% PREC_REC compute precision/recall values
%
%   PREC_REC(PRED,TARGET_POS), where PRED and TARGET_POS are equal-sized 
%   vectors, and TARGET_POS is binary, plots the corresponding precision-recall 
%   graph and the ROC curve.
% 
%   PREC_REC(PRED,TARGET_POS,TARGET_TOTAL) specifies that for all instances 
%   with prediction PRED(i), there are TARGET_POS(i) positive instances and
%   TARGET_TOTAL-TARGET_POS(i) negative ones.
%   
%   PREC_REC(PRED,TARGET_POS,TARGET_TOTAL,NUM_THRESH) additionally specifies
%   the (maximum) number of prediction intervals. Generally, splits are made
%   such that each interval contains about the same number of sample lines.
%   
%   PREC_REC(PRED,TARGET_POS,TARGET_TOTAL,NUM_THRESH,'addToFigure') draws
%   into the current figure, instead of creating a new one.
% 
%   [PREC,TPR,FPR,PRED_THRESH] = PREC_REC(...) does not draw anything, but
%   returns the computed prediction thresholds, along with the respective 
%   precisions, true-positive, and false-positive rates.
%   
%   Example:
%   x1=rand(1000,1);
%   y1=round(x1+0.5*(rand(1000,1)-1));
%   prec_rec(x1,y1);
%   x2=rand(1000,1);
%   y2=round(x2+0.75*(rand(1000,1)-1));
%   prec_rec(x2,y2,[],[],'addToFigure');
%   legend('x1/y1','x2/y2','Location','SouthEast');

%   Stefan Schroedl
%   09/22/2008

% check arguments

if nargin < 2
    error('at least 2 arguments required');
end

[nx,ny]=size(pred);

if (nx~=1 && ny~=1)
    error('first argument must be a vector');
end

[mx,my]=size(pos);
if (mx~=1 && my~=1)
    error('second argument must be a vector');
end
pred = pred(:);
pos = pos(:);

if (length(pos) ~= length(pred))
    error('first two arguments must have same length')
end

if nargin < 3 || isempty(count)
    % set default for total instances
    count = ones(length(pred),1);
else
    [px,py]=size(count);
    if (px~=1 && py~=1)
        error('third argument must be a vector');
    end
    count=count(:);
    if (length(pos) ~= length(count))
        error('first three arguments must have same length')
    end
end

if nargin < 4 || isempty(num_thresh)
    % set default for number of thresholds
    pred_uniq = unique(pred);
    num_thresh = min( floor(length(pred_uniq)/2), 100);
elseif length(num_thresh(:)) ~= 1
    error('fourth argument must be a scalar');
end

keep_figure = 0;
if size(varargin,2)>0 && strcmp(varargin(1),'addToFigure') && ~isempty(get(0,'CurrentFigure'))
        keep_figure = 1;
end
qvals = 0:(1/(num_thresh-1)):1;
pred_thresh = quantile(pred,qvals);
% remove identical bins
pred_thresh = sort(unique(pred_thresh),2,'descend');
total_pos = sum(pos);
total_neg = sum(count-pos);
%score_bins = zeros(length(pred_thresh),1);
prec = zeros(length(pred_thresh),1);
tpr = zeros(length(pred_thresh),1);
fpr = zeros(length(pred_thresh),1);
for i = 1:length(pred_thresh)
    idx = (pred >= pred_thresh(i));
    fpr(i) = sum(count(idx)-pos(idx));
    tpr(i)  = sum(pos(idx))/total_pos;
    prec(i) = sum(pos(idx))/sum(count(idx));
end
fpr = fpr ./total_neg;

if nargout == 0
    
    % draw
    
    if (~keep_figure)
        figure
        subplot(1,2,1);
        hold on
        hold all
        plot(tpr,prec);
        
        xlabel('recall');
        ylabel('precision');
        title('precision-recall graph');
        axis([0 1 0 1]);

        subplot(1,2,2);
        hold on;
        hold all;
        plot(fpr,tpr);
        xlabel('false positive rate');
        ylabel('true positive rate');
        title('roc curve');
        axis([0 1 0 1]);

        % double the width
        rect = get(gcf,'pos');
        rect(3) = 2 * rect(3);
        set(gcf,'pos',rect);
    else
        subplot(1,2,1);
        plot(tpr,prec);

        subplot(1,2,2);
        plot(fpr,tpr);
    end
end