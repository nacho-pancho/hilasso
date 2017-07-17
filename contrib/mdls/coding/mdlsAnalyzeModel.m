% function params=mdlsFitModel(X,D,A)
%
% Given observed samples X, a dictionary and the sparse
% representation A of X in terms of D , this function analyzes the
% data and returns all the relevant model parameters and statistics
%
%
% inputs: 
% X ...... data samples
% D ...... dictionary
% A ...... reconstruction coefficients
%
% outputs:
%
% params . struct with all model parameters
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function params=mdlsAnalyzeModel(X,D,A)
%
% load parameters structure
% 
    A=single(full(A));
    D=single(D);
    X=single(X);
    params=mdlsDefaultParams();
    %
    % empirical approximation error
    %
    E=X-D*A;
    [M,N]=size(X); K=size(D,2);
    E=reshape(E,M*N,1);
    params.muError= mean(E);
    params.varError = var(E);
    [h,c]=hist(E,200);
    E=E-params.muError; % center
    params.histError=struct('freq',h,'center',c);
    lambda=mdlsLaplacianFit(E);
    params.lambdaError=lambda;
    [kappa,betta]=mdlsMOLFit(E);
    params.kappaError=kappa;
    params.betaError=betta;
    clear E;
    %
    % fill in with estimates 
    % for the rec. coefficients
    % first for each row:
    %
    [theta,lambda,mu]=mdlsBerLapFit(A);
    params.muCoef = mu;
    params.lambdaCoef = lambda;
    [theta,kappa,betta]=mdlsBerMOLFit(A);
    params.kappaCoef = kappa;
    params.betaCoef = betta;
    params.thetaCoef = theta;
    Amax=max(max(A)); Amin=min(min(A));
    c=Amin:(Amax-Amin)/100:Amax;
    h=zeros(length(c),K);
    for i=1:K
        h(:,i)=hist(A(i,:),c);
    end
    params.histCoef=struct('freq',h,'center',c);

    scale=1/K;
    NB=25;
    [params.lambdaHist,params.xLambdaHist] = hist(params.lambdaCoef,NB);
    [params.thetaHist,params.xThetaHist] = hist(params.thetaCoef, ...
                                                NB);
    [params.betaHist,params.xBetaHist] = hist(params.betaCoef,NB);
    %
    % fitted Gamma distribution to Laplacian parameter statistics
    %
    [k,b]=mdlsGammaFit(params.lambdaCoef);
    params.kappaCoefGammaPrior = k;
    params.betaCoefGammaPrior = b;
    params.lambdaFit=mdlsGammaDisc(params.xLambdaHist,k,b);
    %params.xLambdaCoef = xLambda;
    %
    % then global ones
    %
    A=reshape(A,1,size(A,1)*size(A,2));
    [theta,lambda]=mdlsBerLapFit(A);
    params.lambdaCoefGlobal = lambda;
    [theta,kappa,betta]=mdlsBerMOLFit(A);
    params.kappaCoefGlobal = kappa;
    params.betaCoefGlobal  = betta;
    params.thetaCoefGlobal = theta;
    [h,c]=hist(A,500);
    params.histCoefGlobal=struct('freq',h,'center',c);
end