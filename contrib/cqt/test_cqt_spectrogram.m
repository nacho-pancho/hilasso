% test constant Q spectrogram code
close all
clear all

[d, fs] = wavread('fornoone2.wav');

% see help cqt_spectrogram for default analysis values
%[CQ fk t] = cqt_spectrogram(d, fs, compute_kernels, Nmax, Nmin, hop, minFreq, maxFreq, bins, scale, window, thresh)

% kernels are computed only once and are stored in a .mat file
compute_kernels = 1; 
tic
[CQ fk t] = cqt_spectrogram(d, fs, compute_kernels);
toc

% then kernels are loaded to save computation
compute_kernels = 0; 
tic
[CQ fk t] = cqt_spectrogram(d, fs, compute_kernels);
toc

figure;
bias = linspace(1,length(fk),length(fk));
bias = bias/bias(end)*10000;
BIAS = repmat(bias',1,length(t));
imagesc(t,fk,20*log10(abs(CQ.*BIAS)+1))
black = abs(1-gray);
colormap(black);set(gca,'YDir','normal');
title('Constant Q transform');xlabel('Time (s)');ylabel('Frequency (Hz)');

