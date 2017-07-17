 --------------------------------------------------------------------------
  PACKAGE: glmnet
 --------------------------------------------------------------------------

  TITLE: Lasso and elastic-net regularized generalized linear models
  
  DESCRIPTION:
    This package fits lasso and elastic-net model paths for regression, logistic and multinomial regression
    using coordinate descent. The algorithm is extremely fast, and exploits sparsity in the input x
    matrix where it exists. A variety of predictions can be made from the fitted models.

  LICENSE: GPL-2
 
  DATE: 14 Jul 2009

  FUNCTIONS:  
     glmnet
     glmnetSet
     glmnetPlot
     glmnetPredict
     glmnetPrint
     glmnetCoef
 
  AUTHORS:
     Algorithm was designed by Jerome Friedman, Trevor Hastie and Rob Tibshirani 
     Fortran code was written by Jerome Friedman 
     R wrapper (from which the MATLAB wrapper was adapted) was written by Trevor Hasite
     MATLAB wrapper was written and maintained by Hui Jiang, jiangh@stanford.edu 
     Department of Statistics, Stanford University, Stanford, California, USA.
 
  REFERENCES:
     Friedman, J., Hastie, T. and Tibshirani, R. (2009)
     Regularization Paths for Generalized Linear Models via Coordinate Descent.
     To appear, Journal of Statistical Software
 
  URL: http://www-stat.stanford.edu/~hastie/Papers/glmnet.pdf

  EXAMPLES:
     x=randn(100,20);
     y=randn(100,1);
     g2=randsample(2,100,true);
     g4=randsample(4,100,true);
     fit1=glmnet(x,y);
     glmnetPredict(fit1, 'link', x(1:5,:), [0.01,0.005]')
     glmnetPredict(fit1, 'coefficients');
     glmnetPlot(fit1, 'lambda');
     fit2=glmnet(x,g2,'binomial');
     glmnetPredict(fit2,'response',x(2:5,:))
     glmnetPredict(fit2,'nonzero')
     fit3=glmnet(x,g4,'multinomial');
     glmnetPredict(fit3,'response',x(1:3,:),0.01)

 --------------------------------------------------------------------------
  glmnet.m: fit an elasticnet model path
 --------------------------------------------------------------------------
 
  DESCRIPTION:
     Fit a regularization path for the elasticnet at a grid of values for
     the regularization parameter lambda. Can deal with all shapes of data.
     Fits linear, logistic and multinomial regression models.
 
  USAGE: 
     fit = glmnet(x, y)
     fit = glmnet(x, y, family, options)
 
  EXTERNAL FUNCTIONS:
  options         = glmnetSet;                  provided with glmnet.m
 
  INPUT ARGUMENTS:
  x           Input matrix, of dimension nobs x nvars; each row is an
              observation vector. Can be in sparse column format.
  y           Response variable. Quantitative for family =
              'gaussian'. For family = 'binomial' should be a two-column
              matrix of counts or proportions. For family = 'multinomial',
              should be a matrix with nc columns of counts or proportions.
  family      Reponse type. (See above). Default is 'gaussian'. 
  options     A structure that may be set and altered by glmnetSet (type 
              help glmnetSet).
 
  OUTPUT ARGUMENTS:
  fit         A structure.
  fit.a0      Intercept sequence of length length(fit.lambda). 
  fit.beta    For "elnet" and "lognet" models, a nvars x length(lambda) 
              matrix of coefficients. For "multnet", a list of nc such
              matrices, one for each class.
  fit.lambda  The actual sequence of lambda values used.
  fit.dev     The fraction of (null) deviance explained (for "elnet", this 
              is the R-square).
  fit.nulldev Null deviance (per observation).
  fit.df      The number of nonzero coefficients for each value of lambda.
              For "multnet", this is the number of variables with a nonzero
              coefficient for any class.
  fit.dfmat   For "multnet" only. A matrix consisting of the number of
              nonzero coefficients per class.
  fit.dim     Dimension of coefficient matrix (ices).
  fit.npasses Total passes over the data summed over all lambda values.
  fit.jerr    Error flag, for warnings and errors (largely for internal 
              debugging).
  fit.type    Type of regression - internal usage.
 
  DETAILS:
     The sequence of models implied by lambda is fit by coordinate descent.
     For family = 'gaussian' this is the lasso sequence if alpha = 1, else
     it is the elasticnet sequence. For family = 'binomial' or family =
     "multinomial", this is a lasso or elasticnet regularization path for
     fitting the linear logistic or multinomial logistic regression paths.
     Sometimes the sequence is truncated before options.nlambda values of
     lambda have been used, because of instabilities in the logistic or
     multinomial models near a saturated fit. glmnet(..., family =
     'binomial') fits a traditional logistic regression model for the
     log-odds. glmnet(..., family = 'multinomial') fits a symmetric
     multinomial model, where each class is represented by a linear model
     (on the log-scale). The penalties take care of redundancies. A
     two-class "multinomial" model will produce the same fit as the
     corresponding "binomial" model, except the pair of coefficient
     matrices will be equal in magnitude and opposite in sign, and half the
     "binomial" values. Note that the objective function for
     "gaussian" is 
                 1 / (2 * nobs) RSS + lambda * penalty
     , and for the logistic models it is 
                 1 / nobs - loglik + lambda * penalty
 
  LICENSE: GPL-2
 
  DATE: 14 Jul 2009
 
  AUTHORS:
     Algorithm was designed by Jerome Friedman, Trevor Hastie and Rob Tibshirani 
     Fortran code was written by Jerome Friedman 
     R wrapper (from which the MATLAB wrapper was adapted) was written by Trevor Hasite
     MATLAB wrapper was written and maintained by Hui Jiang, jiangh@stanford.edu 
     Department of Statistics, Stanford University, Stanford, California, USA.
 
  REFERENCES:
     Friedman, J., Hastie, T. and Tibshirani, R. (2009)
     Regularization Paths for Generalized Linear Models via Coordinate Descent.
     To appear, Journal of Statistical Software
 
  SEE ALSO:
     glmnetSet, glmnetPrint, glmnetPlot, glmnetPredict and glmnetCoef methods.
  
  EXAMPLES:
     x=randn(100,20);
     y=randn(100,1);
     g2=randsample(2,100,true);
     g4=randsample(4,100,true);
     fit1=glmnet(x,y);
     glmnetPrint(fit1);
     glmnetCoef(fit1,0.01) % extract coefficients at a single value of lambda
     glmnetPredict(fit1,'response',x(1:10,:),[0.01,0.005]') % make predictions
     fit2=glmnet(x,g2,'binomial');
     fit3=glmnet(x,g4,'multinomial');
 
  DEVELOPMENT: 14 Jul 2009: Original version of glmnet.m written.
  
 --------------------------------------------------------------------------
  glmnetPlot.m: plot coefficients from a "glmnet" object
 --------------------------------------------------------------------------
 
  DESCRIPTION:
     Produces a coefficient profile plot fo the coefficient paths for a
     fitted "glmnet" object.
 
  USAGE: 
     glmnetPlot(fit);
     glmnetPlot(fit, xvar);
     glmnetPlot(fit, xvar, label);
 
  INPUT ARGUMENTS:
  x           fitted "glmnet" model.
  xvar        What is on the X-axis. "norm" plots against the L1-norm of
              the coefficients, "lambda" against the log-lambda sequence,
              and "dev" against the percent deviance explained.
  label       if TRUE, label the curves with variable sequence numbers.
 
  DETAILS:
     A coefficient profile plot is produced. If x is a multinomial model, a
     coefficient plot is produced for each class.
 
  LICENSE: GPL-2
 
  DATE: 14 Jul 2009
 
  AUTHORS:
     Algorithm was designed by Jerome Friedman, Trevor Hastie and Rob Tibshirani 
     Fortran code was written by Jerome Friedman 
     R wrapper (from which the MATLAB wrapper was adapted) was written by Trevor Hasite
     MATLAB wrapper was written and maintained by Hui Jiang, jiangh@stanford.edu 
     Department of Statistics, Stanford University, Stanford, California, USA.
 
  REFERENCES:
     Friedman, J., Hastie, T. and Tibshirani, R. (2009)
     Regularization Paths for Generalized Linear Models via Coordinate Descent.
     To appear, Journal of Statistical Software
 
  SEE ALSO:
     glmnet, glmnetSet, glmnetPrint, glmnetPredict and glmnetCoef methods.
  
  EXAMPLES:
     x=randn(100,20);
     y=randn(100,1);
     g2=randsample(2,100,true);
     g4=randsample(4,100,true);
     fit1=glmnet(x,y);
     glmnetPlot(fit1);
     glmnetPlot(fit1, 'lambda', true);
     fit3=glmnet(x,g4,'multinomial');
     glmnetPlot(fit3);
 
  DEVELOPMENT: 14 Jul 2009: Original version of glmnet.m written.

 --------------------------------------------------------------------------
  glmnetPredict.m: make predictions from a "glmnet" object.
 --------------------------------------------------------------------------
 
  DESCRIPTION:
     Similar to other predict methods, this functions predicts fitted
     values, logits, coefficients and more from a fitted "glmnet" object.
 
  USAGE: 
     glmnetPredict(object)
     glmnetPredict(object, type)
     glmnetPredict(object, type, newx)
     glmnetPredict(object, type, newx, s)
 
  INPUT ARGUMENTS:
  fit         Fitted "glmnet" model object.
  type        Type of prediction required. Type "link" gives the linear
              predictors for "binomial" or "multinomial" models; for
              "gaussian" models it gives the fitted values. Type "response"
              gives the fitted probabilities for "binomial" or
              "multinomial"; for "gaussian" type "response" is equivalent
              to type "link". Type "coefficients" computes the coefficients
              at the requested values for s. Note that for "binomial"
              models, results are returned only for the class corresponding
              to the second level of the factor response. Type "class"
              applies only to "binomial" or "multinomial" models, and
              produces the class label corresponding to the maximum
              probability. Type "nonzero" returns a list of the indices of
              the nonzero coefficients for each value of s.
  newx        Matrix of new values for x at which predictions are to be 
              made. Must be a matrix; This argument is not used for 
              type=c("coefficients","nonzero")
  s           Value(s) of the penalty parameter lambda at which predictions
              are required. Default is the entire sequence used to create
              the model.
 
  DETAILS:
     The shape of the objects returned are different for "multinomial"
     objects. glmnetCoef(fit, ...) is equivalent to glmnetPredict(fit, "coefficients", ...)
 
  LICENSE: GPL-2
 
  DATE: 14 Jul 2009
 
  AUTHORS:
     Algorithm was designed by Jerome Friedman, Trevor Hastie and Rob Tibshirani 
     Fortran code was written by Jerome Friedman 
     R wrapper (from which the MATLAB wrapper was adapted) was written by Trevor Hasite
     MATLAB wrapper was written and maintained by Hui Jiang, jiangh@stanford.edu 
     Department of Statistics, Stanford University, Stanford, California, USA.
 
  REFERENCES:
     Friedman, J., Hastie, T. and Tibshirani, R. (2009)
     Regularization Paths for Generalized Linear Models via Coordinate Descent.
     To appear, Journal of Statistical Software
 
  SEE ALSO:
     glmnet, glmnetSet, glmnetPrint, glmnetPlot and glmnetCoef methods.
  
  EXAMPLES:
     x=randn(100,20);
     y=randn(100,1);
     g2=randsample(2,100,true);
     g4=randsample(4,100,true);
     fit1=glmnet(x,y);
     glmnetPredict(fit1,'link',x(1:5,:),[0.01,0.005]') % make predictions
     glmnetPredict(fit1,'coefficients')
     fit2=glmnet(x,g2,'binomial');
     glmnetPredict(fit2, 'response', x(2:5,:))
     glmnetPredict(fit2, 'nonzero')
     fit3=glmnet(x,g4,'multinomial');
     glmnetPredict(fit3, 'response', x(1:3,:), 0.01)
 
  DEVELOPMENT: 14 Jul 2009: Original version of glmnet.m written.
  
 --------------------------------------------------------------------------
  glmnet.m: print a glmnet object
 --------------------------------------------------------------------------
 
  DESCRIPTION:
     Print a summary of the glmnet path at each step along the path.
 
  USAGE: 
     glmnetPrint(fit)
 
  INPUT ARGUMENTS:
  fit         fitted glmnet object
 
  DETAILS:
     Three-column matrix with columns Df, dev and Lambda is printed. The Df
     column is the number of nonzero coefficients (Df is a reasonable name
     only for lasso fits). dev is the percent deviance explained (relative
     to the null deviance).
 
  LICENSE: GPL-2
 
  DATE: 14 Jul 2009
 
  AUTHORS:
     Algorithm was designed by Jerome Friedman, Trevor Hastie and Rob Tibshirani 
     Fortran code was written by Jerome Friedman 
     R wrapper (from which the MATLAB wrapper was adapted) was written by Trevor Hasite
     MATLAB wrapper was written and maintained by Hui Jiang, jiangh@stanford.edu 
     Department of Statistics, Stanford University, Stanford, California, USA.
 
  REFERENCES:
     Friedman, J., Hastie, T. and Tibshirani, R. (2009)
     Regularization Paths for Generalized Linear Models via Coordinate Descent.
     To appear, Journal of Statistical Software
 
  SEE ALSO:
     glmnet, glmnetSet, glmnetPlot, glmnetPredict and glmnetCoef methods.
  
  EXAMPLES:
     x=randn(100,20);
     y=randn(100,1);
     fit1=glmnet(x,y);
     glmnetPrint(fit1);
 
  DEVELOPMENT: 14 Jul 2009: Original version of glmnet.m written.
  
 --------------------------------------------------------------------------
  glmnetCoef.m: print coefficients from a "glmnet" object
 --------------------------------------------------------------------------
 
  USAGE: 
     glmnetCoef(fit);
     glmnetCoef(fit, s);
 
  DETAILS:
     glmnetCoef(fit, s) is equivalent to glmnetPredict(fit, "coefficients", [], s)
     See glmnetPredict for more details.
 
  LICENSE: GPL-2
 
  DATE: 14 Jul 2009
 
  AUTHORS:
     Algorithm was designed by Jerome Friedman, Trevor Hastie and Rob Tibshirani 
     Fortran code was written by Jerome Friedman 
     R wrapper (from which the MATLAB wrapper was adapted) was written by Trevor Hasite
     MATLAB wrapper was written and maintained by Hui Jiang, jiangh@stanford.edu 
     Department of Statistics, Stanford University, Stanford, California, USA.
 
  REFERENCES:
     Friedman, J., Hastie, T. and Tibshirani, R. (2009)
     Regularization Paths for Generalized Linear Models via Coordinate Descent.
     To appear, Journal of Statistical Software
 
  SEE ALSO:
     glmnet, glmnetSet, glmnetPrint, glmnetPredict and glmnetPlot methods.
  
  EXAMPLES:
     x=randn(100,20);
     y=randn(100,1);
     fit1=glmnet(x,y);
     glmnetCoef(fit1,0.01) % extract coefficients at a single value of lambda
 
  DEVELOPMENT: 14 Jul 2009: Original version of glmnet.m written.
  
  