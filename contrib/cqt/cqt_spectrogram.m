function [CQ fk t] = cqt_spectrogram(d, fs, compute_kernels, Nmax, Nmin, hop, minFreq, maxFreq, bins, scale, window, thresh)
% function [CQ fk t] = cqt_spectrogram(d, fs, compute_kernels, Nmax, Nmin, hop, minFreq, maxFreq, bins, scale, window, thresh)
%
% Computes a constant Q spectrogram using the algorithm described in:
%   "An efficient algorithm for the calculation of a constant Q transform"
%   Judith C. Brown and Miller S. Puckette, Journal of the Acoustic Society of America, 1992
%
% This code is based on an implementation by Benjamin Blankertz but with
% some corrections and changes:
% "The Constant Q Transform", Benjamin Blankertz
% http://wwwmath.uni-muenster.de/logik/Personen/blankertz/constQ/constQ.html},
% 
% input:
%               d - audio samples
%              fs - sampling frequency (default values assume fs = 44100)
% compute_kernels - whether to compute the kernles or load them from a
%                   previously created file (0/1 - default 1)
%            Nmax - max number of samples in the frame (default 4096)
%            Nmin - min number of samples in the frame (default  512)
%             hop - hop size in samples (default 256)
%         minFreq - min frequency bin to compute (default   60)
%         maxFreq - min frequency bin to compute (default 6000)
%            bins - number of bins in an octave (12 = semi-tone, 24 =
%                   quarter-tone, etc - default 48)
%           scale - type of bin frequency distribution (linear or
%                   logarithmic - default linear)
%          window - type of window used (hamming or hann - default hamming)
% 
% output:
%              CQ - spectrogram computed using the constant-Q transform 
%                   each column contains an estimate of the short-term
%                   spectrum of d, time increases across the columns and 
%                   frequency increases down the rows
%              fk - vector of frequencies in Hz 
%               t - vector of time instants in seconds, each value
%                   corresponds to the center of each segment

if nargin< 2;             fs   =     44100; end
if nargin< 3;  compute_kernels =         1; end
if nargin< 4;             Nmax =      4096; end
if nargin< 5;             Nmin =       512; end
if nargin< 6;             hop  =       256; end
if nargin< 7;          minFreq =        60; end
if nargin< 8;          maxFreq =      6000; end
if nargin< 9;             bins =        48; end
if nargin<10;            scale =  'linear'; end
if nargin<11;           window = 'hamming'; end
if nargin<12;           thresh =    0.0054; end  % for Hamming window

% compute or load kernels
if compute_kernels
    disp('Computing kernels ...')
    [sparKernel fk] = sparseKernelComputation(minFreq, maxFreq, bins, fs, Nmax, Nmin, scale, window, thresh);
    save kernels sparKernel fk;
else
    disp('Loading kernels ...')
    load kernels
end

% compute constant-Q spectrogram
K = length(fk);
len = 2^nextpow2(Nmax);
M = floor(length(d)/hop);
d = [d; zeros(len-mod(length(d),hop),1)];
CQ = zeros(K,M);
S = size(sparKernel,1);
d = d(:);

for i=1:M
    tini = (i-1)*hop+1;
    x = d(tini:tini+len-1)';
    cq = fft(x,S) * sparKernel;
    CQ(:,i) = cq';
end

t = (hop:hop:M*hop)/fs;

end

function [sparKernel fk specKernels tempKernels lens Qs] = sparseKernelComputation(minFreq, maxFreq, bins, fs, Nmax, Nmin, scale, window, thresh)

% function [sparKernel fk specKernels tempKernels lens Qs] = sparseKernelComputation(minFreq, maxFreq, bins, fs, Nmax, Nmin, scale, window, thresh)
%
% Function to compute sparse kernels for the efficient implementation of a
% constant-q transform. Frequency kernels are near zero over most of
% the spectrum. For this reason a sparse representation is used in order to
% save computations. 
% 
if nargin<7;  scale = 'linear'; end
if nargin<8; window = 'hamming'; end
if nargin<9; thresh = 0.0054; end   % for Hamming window (may be different for Hann)

% Q factor
Q = 1/(2^(1/bins)-1);

% type of frequency scale: linear or logarithmic
if strcmp(scale, 'linear')
    fk = minFreq:fs/Nmax:maxFreq;
    K = length(fk);
else
    K = ceil( bins * log2(maxFreq/minFreq) );
    fk = minFreq*2.^(((1:1:K)-1)/bins);
end

% window lengths (forced to be even, and bounded between Nmax and Nmin) 
lens = ceil( Q * fs ./ fk);
lens(mod(lens,2)~=0) = lens(mod(lens,2)~=0) + 1;
lens(lens > Nmax) = Nmax; % Nmax and Nmin are assumed to be even
lens(lens < Nmin) = Nmin;

% effective quality factors (forced to be integer)
Qs = fk.*lens/fs;
                       
sparKernel = [];
if nargout > 2
    tempKernels = zeros(Nmax, K);
    specKernels = zeros(Nmax, K);
end

for k= K:-1:1;
    len = lens(k);
    Q   = Qs(k);
    % window type: hamming or hann
    if strcmp(window, 'hamming')
        win = hamming(len,'symmetric')/len;
    else
        win = hann(len,'symmetric')/len;
    end
    % build temporal kernel
    kern  = exp(2*pi*i*Q*(0:len)'/(len));
    kern = kern(1:end-1);
    midp = floor(len/2);
    tempKernel =  fftshift([win(midp+1:end) .* kern(midp+1:end); zeros(Nmax-len,1) ; win(1:midp) .* kern(1:midp)]);  
    % spectral kernel
    specKernel = fft(tempKernel);
    % thresholding 
    specKernel(abs(specKernel)<=thresh) = 0;
    sparKernel = sparse([specKernel sparKernel]);
    % build sparse kernel
    if nargout > 2
        tempKernels(:,k) = tempKernel;
        specKernels(:,k) = specKernel;
    end
end

% spectral kernels
sparKernel = conj(sparKernel) / Nmax;

end
