x2=textread('x2.txt');
g4=textread('g4.txt');
fit3=glmnet(x2,g4,'multinomial');
glmnetPredict(fit3,'response',x2(2:5,:))