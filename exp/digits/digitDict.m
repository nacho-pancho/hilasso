
load data/usps/usps_train.mat

inipar = mdlsDictIni();
inipar.X = data.X;
inipar.labels = data.label;
inipar.method = 'patches';
inipar.sigma = 0.01;
inipar.K = 300;
D0     = mdlsDictIni(inipar);
clear inipar;

learnpar = mdlsLearnModel();
learnpar.base_name = 'usps';
learnpar.output_dir= 'results/dictionaries/usps';
system('mkdir -p results/dictionaries/usps');
learnpar.training_data = data.X;
learnpar.training_labels = data.label;
learnpar.D0 = D0;

D = mdlsLearnModel(learnpar);
clear learnpar
