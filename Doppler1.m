clear all, close all, clc
SatPassFile = uigetfile('*.wav')

FileCode = textscan(SatPassFile,'%s%s%s%s%s','delimiter','_');

% Time of recording
date = FileCode{2}{1};
time = FileCode{3}{1}(1:end-1);

% Central Frequency
F_central = FileCode{4}{1};



[y,FS,NBITS,OPTS]=wavread(SatPassFile);

signal = complex(y(:,1),y(:,2)); 
% first channel, what is the second? stereo...? complex sampling?
%Stereo most likely!

TS = 1/FS; % Sampling period

% Split the signal in blocks of LT seconds
LT = 0.005; 
L = ceil(1*FS); % Number of samples to analyse
NFFT = 2^nextpow2(L);


% Reduce spectral resolution by factor DS
DS = 100;
SPEC = zeros(floor(NFFT/DS),floor(numel(signal)/L));

% Compute spectra for the time windows defined by LT
for kLT = 1 : floor(numel(signal)/L);
    
IndSignal = ((L * (kLT-1)) : (L * kLT)) + 1;

S = fft(signal(IndSignal),NFFT)/L;
S = [S(NFFT/2+2:end);S(1:NFFT/2+1)];
f = FS/2*linspace(-1,1,NFFT);


S_lowres = decimate(abs(S),DS);
f_lowres = decimate(f,DS);

SPEC(1:numel(S_lowres),kLT) = S_lowres; 

plot(f,abs(S))
drawnow
end
t = TS * (0:(kLT-1));

% Noise estimate and SNR
Noise = mean(mean(SPEC(1:30,:)));
SNR = 10*log10(SPEC/Noise);
SNR(imag(SNR)~=0) = 0;


% Stackplot of spectra
surf(t,f_lowres*1e-3,SNR)
shading flat
view(2)
cb = colorbar;
set(get(cb,'ylabel'),'String', 'SNR [dB]');
caxis([0,max(caxis)])
title([date,' ',time,' Central Freq. ',F_central])
xlabel('Seconds [s]')
ylabel('Doppler [kHz]')
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of spectrum
%Start of Doppler measurements
%SNR is decibel of signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doppler1=zeros(length(SNR),length(t));
SNRsizearray=size(SNR);
SNRsize=SNRsizearray(1)*SNRsizearray(2);
fdoppler=zeros(length(SNR),length(t));
for j=1:length(t);
    fdoppler(:,j)=f_lowres';
end

%Matrix with Frequency corresponding to decibel higher than 10 and column
%corresponding to the time index.
for i=1:SNRsize;
    if SNR(i)>10
       doppler1(i)=fdoppler(i)*1e-3;
    end
end


%NEXT How to plot several frequencies corresponding to times?
