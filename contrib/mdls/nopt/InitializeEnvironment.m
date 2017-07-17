cd ..
global NOPT;
NOPT=pwd;
cd code
addpath('./experiments',          '-end');
%addpath('../lib/twist',      '-end');
addpath('../lib/julien',      '-end');
addpath('../lib/ksvd',        '-end');
addpath('../dictionaries',    '-end');
addpath('../images',          '-end');
addpath('./util',             '-end');
addpath('./priors',           '-end');
addpath('./sandbox',          '-end');
