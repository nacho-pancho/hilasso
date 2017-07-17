% function stats = mdlsCoeffStatsLearn(D,X,e[,stats])
%
% This function takes a dictionary D and data X, and optionally already
% computed statistics from a previous run, stats, and constructs 
% a probability model for sparse recovery coefficients based on those statistics.
% The reconstruction of the data X is exact up to a maximum distortion e
% (squared error).
% The data is assumed normalized, that is, each row must have zero mean
% and norm 1.
% 
% INPUT
%   D ............. (MxK) dictionary
%   X ............. (MxN) training data, *each column assumed to be zero mean and norm 1*
%   e ............. (1x1) maximum allowable squared distortion. If empty
%                         or set to <= zero it will take a default value of 1e-6.
%   stats ......... (struct) statistics learned from a previous run or
%                            during mdlsDictLearn. Fields:
%       AAt ............... empirical covariance matrix of coefficients
%       XAt ............... empirical cross-correlation of data and 
%                             coeffs.
%       rho ............... frequency of usage of each atom.
%       xrho .............. co-occurence matrix of atoms, that is, 
%                             abs(sgn(A))*abs(sgn(A))^T
%       N ................. effective number of samples processed.

% OUTPUT
%
% stats ........... Updated statistics.
%
function stats = mdlsCoeffStatsLearn(D,X,e,stats)

%======================================================================
% INITIALIZATION
%======================================================================
[M,K]=size(D);
if ~exist('stats','var')
    stats = struct();
    stats.AAt = zeros(K,K);
    stats.XAt = zeros(M,K);
    stats.rho = zeros(K,1);
    stats.xrho= zeros(K,K);
    stats.N   = 0;
end
if ~exist('e','var')
    e = 0;
end
if ~exist('X','var')
    X =[];
end
if isempty(X)
    return;
end
[Mx,N]= size(X);
if M ~= Mx
  error('Incompatible dictionary for data');
end
if e <= 0
    e=1e-6;
end
%
%======================================================================
% L1 REGULARIZED CODING
%======================================================================
%
s = min(K,M-1);
l1_mode = 1;
A = scLasso(X,D,s,e,l1_mode);
%
%======================================================================
% UPDATE STATISTICS
%======================================================================
%
stats.AAt = stats.AAt + A*A';
stats.XAt = stats.XAt + X*A';
A = abs(sign(A));
stats.rho = stats.rho + sum(A,2);
stats.xrho = stats.xrho + A*A';
stats.N = stats.N + N;
