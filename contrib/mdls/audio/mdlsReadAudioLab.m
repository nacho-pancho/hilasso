%
% function L=mdlsReadAudioLab(fname,Fs)
%
% this will read a .lab file with annotated audio
% segments and create an 1xN matrix where each 
% sample is labeled according to the lab file.
% the sample frequency, which is needed, is given in Fs
% in hertz.
%
function L=mdlsReadAudioLab(fname,Fs)
    fid=fopen(fname,'r');
    A = textscan(fid,'%n%n%s');
    fclose(fid);
    N=ceil(Fs*A{2}(end)); % total number of samples
                          %fprintf('Total number of samples: %d\n',N);
    j=2;
    NR = length(A{1});
    si = 1;
    L=zeros(1,N);
    for i=1:NR
        label = A{3}{i};
        sf = floor( A{2}(i) * Fs );
        L(si:sf)=label(1);
        si = sf+1;
    end 
    L=L(1:sf);
end