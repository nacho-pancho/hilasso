x = [1 1; 2 2; 3 3];
y = [1 3 2]';
fit1 = glmnet(x, y);
glmnetPrint(fit1);
glmnetPlot(fit1);
g2 = [1 1 2]';
g3 = [1 2 3]';
fit2 = glmnet(x, g2, 'binomial');
fit3 = glmnet(x, g3, 'multinomial');
glmnetPrint(fit2);
glmnetPrint(fit3);
glmnetPredict(fit1,'coefficients')
glmnetPredict(fit1,'link',x(1:2,:),[0.01, 0.005]')
glmnetCoef(fit1, 0.01)
glmnetPredict(fit2,'response',x(1:2,:))
glmnetPredict(fit2,'nonzero')
glmnetPredict(fit3,'response',x(1:3,:),0.01)

x = randn(100,20);
y = randn(100,1);
fit1 = glmnet(x, y);
glmnetPrint(fit1);
glmnetPlot(fit1);
g2 = randsample(2,100,true);
g4 = randsample(4,100,true);
fit2 = glmnet(x, g2, 'binomial');
fit3 = glmnet(x, g4, 'multinomial');
glmnetPrint(fit2);
glmnetPrint(fit3);
glmnetPredict(fit1,'coefficients')
glmnetPredict(fit1,'link',x(2:5,:),[0.01, 0.005]')
glmnetCoef(fit1, 0.01)
glmnetPredict(fit2,'response',x(2:5,:))
glmnetPredict(fit2,'nonzero')
glmnetPredict(fit3,'response',x(1:3,:),0.01)
glmnetPlot(fit2);
glmnetPlot(fit3);