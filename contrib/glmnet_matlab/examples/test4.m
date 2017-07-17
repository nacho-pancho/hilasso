x1=textread('x1.txt');
g2=textread('g2.txt');
fit2=glmnet(x1,g2,'binomial');
glmnetPredict(fit2,'response',x1(2:5,:))