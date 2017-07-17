%
% function mdlsPrintModel(modelParams)
%
% Name: mdlsPrintModel
%
% Category: utility
%
% Description: print description of a sparse model
%
% Input:
% modelParams ... struct with model parameters. See mdlsDefaultParams
%
%
% Output:
% none
%
% Author: Ignacio Ramirez <ignacio.ramirez at gmail dot com>
%
% Version: $Id$
%
function mdlsPrintModel(m,name)
    format short;
    fprintf('MODEL:%s\n', name);
    fprintf('A: Ber:theta=%5.2f\tLap:lambda=%5.2f',m.thetaCoefGlobal,m.lambdaCoefGlobal);
    fprintf('\tMOL:kappa=%5.2f,beta=%5.2f\n', ...
            m.kappaCoefGlobal,m.betaCoefGlobal);
    fprintf('lambda: kappa=%5.2f,beta=%5.2f. chi-square score=%f\n', ...
            m.kappaCoefGammaPrior,m.betaCoefGammaPrior,...
            sum( (m.lambdaHist-m.lambdaFit).^2./m.lambdaFit ));
    fprintf('E: Gauss:var=%8.6f\tLap:lambda=%5.2f',m.varError,m.lambdaError);
    fprintf('\tMOL:kappa=%5.2f,beta=%5.2f\n',m.kappaError,m.betaError);
end