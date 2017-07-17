%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
%
global WORKSPACE DATA IMAGES RESULTS;
WORKSPACE=[pwd() '/'];
DATA=[WORKSPACE 'data/'];
IMAGES=[DATA 'images/'];
RESULTS=[WORKSPACE 'results/'];
% 
% my library
%
addpath([WORKSPACE 'contrib/mdls'],'-end');
addpath([WORKSPACE 'contrib/mdls/audio'],'-end'); % audio-related functions
addpath([WORKSPACE 'contrib/mdls/image'],'-end'); % image-related functions
addpath([WORKSPACE 'contrib/mdls/util'],'-end');  % misc utilities
addpath([WORKSPACE 'contrib/mdls/coding'],'-end');% sparse coding 
addpath([WORKSPACE 'contrib/mdls/modeling'],'-end');% sparse model learning
addpath([WORKSPACE 'contrib/mdls/priors'],'-end'); % custom probability models
addpath([WORKSPACE 'contrib/mdls/infotheo'],'-end'); % information-theoretic
addpath([WORKSPACE 'contrib/mdls/mdl'],'-end');      % MDL-related
addpath([WORKSPACE 'contrib/mdls/classification'],'-end'); 
addpath([WORKSPACE 'contrib/mdls/clustering'],'-end');
addpath([WORKSPACE 'contrib/mdls/spams_wrapper'],'-end');
%addpath([WORKSPACE 'contrib/mdls/hilasso'],'-end');
%addpath([WORKSPACE 'contrib/mdls/hilasso/SpaRSA_2'],'-end');
addpath([WORKSPACE 'src/util'],'-end');
addpath([WORKSPACE 'src/nacho'],'-end');
addpath([WORKSPACE 'src/nacho/SpaRSA_2'],'-end');
%
% contributed code
%
addpath([WORKSPACE 'contrib'],'-end');
%
% SIFT
%
addpath( [WORKSPACE 'contrib/ScSPM'],'-end');   % 
%
% Precision-recall curves
%
addpath( [WORKSPACE 'contrib/prec_rec'],'-end');   % 
%
% Graphcut/spectral clustering
%
addpath( [WORKSPACE 'contrib/graph'],'-end');   % 
addpath( [WORKSPACE 'contrib/speclu'],'-end'); % Wen-Yen Chen
%
% Modified Cholesky implementation
%
addpath([WORKSPACE 'contrib/mdoc'],'-end');  
%
% Sparse Modeling library SPAMS
%
if isequal(mexext(),'mexglx')
  addpath([WORKSPACE 'contrib/SPAMS/release/mkl32'],'-end');
elseif isequal(mexext(),'mexa64')
  addpath([WORKSPACE 'contrib/SPAMS/release/mkl64'],'-end');
elseif isequal(mexext(),'mexmaci')
  addpath([WORKSPACE 'contrib/SPAMS/release/macmkl32'],'-end');
else
  warning('No implementation available for SPAMS');
end
%
% constant quality transform
%
addpath([WORKSPACE 'contrib/cqt'],'-end');   
%
% CVX
%
addpath( [WORKSPACE 'contrib/cvx/builtins'],'-end');
addpath( [WORKSPACE 'contrib/cvx/commands'],'-end');
addpath ([WORKSPACE 'contrib/cvx/functions'],'-end');
addpath ([WORKSPACE 'contrib/cvx/lib'],'-end');
addpath ([WORKSPACE 'contrib/cvx/structures'],'-end');
addpath ([WORKSPACE 'contrib/cvx'],'-end');
%
% VOICEBOX - audio analysis library
%
%addpath( [WORKSPACE 'contrib/voicebox'],'-end');   % 
%
% scripts for experiments
%
[status,expdir1]=system(['find exp -type d | grep -v .svn | grep -v /old/']);
expdirs = mdlsSplitCell(expdir1,char(10));
for i=1:length(expdirs)
    addpath([WORKSPACE expdirs{i}],'-end');
end
redblue=zeros(64,3);
%
% colormaps for seeing negative/positive structures
% (dictioanries)
%
grad=(1:32)/32;
redblue(32:-1:1,1)=grad;
redblue(33:64,3)=grad;

redgreen=zeros(64,3);
redgreen(32:-1:1,1)=grad;
redgreen(33:64,2)=grad;

redcyan=redgreen;
redcyan(33:64,2)=grad;

%
% homework scripts
%

