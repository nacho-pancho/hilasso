x=textread('x.txt');
y=textread('y.txt');
fit1=glmnet(x,y);
glmnetPrint(fit1);
glmnetPlot(fit1);

b1=glmnetCoef(fit1);

f1=glmnetPredict(fit1,'link',x);
