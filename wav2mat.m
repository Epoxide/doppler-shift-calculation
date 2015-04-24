clear all, close all, clc
[SatPassFile,Pathname] = uigetfile('*.wav');

FileCode = textscan(SatPassFile,'%s%s%s%s%s','delimiter','_');

% Time of recording
date = FileCode{2}{1};
time = FileCode{3}{1}(1:end-1);

starttimemes = datenum([str2num(date(1:4)), str2num(date(5:6)),...
    str2num(date(7:8)), str2num(time(1:2)), str2num(time(3:4)),...
    str2num(time(5:6))]);


% Central Frequency
F_central = str2num(FileCode{4}{1}(1:end-3))*1e3;


% Get Data Info
FileInfo = audioinfo([Pathname,SatPassFile]);
FS = FileInfo.SampleRate;

TS = 1/FS; % Sampling period

% Split the signal in blocks of LT seconds
LT = 0.5;
L = ceil(1*FS*LT); % Number of samples to analyse
NFFT = 2^nextpow2(L);

% Reduce spectral resolution by factor DS
DS = 10;

% Allocate memory for 10 min of data
SPEC = zeros(floor(NFFT/DS),floor(FileInfo.TotalSamples/L));
windowf = hamming(L);

% Compute spectra for the time windows defined by LT
for kLT = 1 : floor(FileInfo.TotalSamples/L)
    % Read the data from file
    IndSignal = (L * (kLT-1) + 1 ) : ((L * kLT));
    
    [y,FS]=audioread([Pathname,SatPassFile],[IndSignal(1) IndSignal(end)]);
    signal = complex(y(:,1),y(:,2));
    % complex sampling?
    
    
    S = abs(real(fft(windowf.*signal,NFFT)/L));
    S = flipud([S(NFFT/2+2:end);S(1:NFFT/2+1)]);
    f = FS/2*linspace(-1,1,NFFT) + F_central;
    
    
    S_lowres = decimate(S,DS);
    f_lowres = mean(reshape(f(1:floor(numel(f)/DS)*DS),DS,floor(numel(f)/DS)),1);
    %decimate(f*1e-6,DS)*1e6;
    
    SPEC(1:numel(S_lowres),kLT) = S_lowres;
    
end

% Time in seconds from start
t =  TS * (0:(kLT-1))*L;


% Remove Background
BKGR = mean(SPEC(:,1:2),2)*1;
SPEC = SPEC - repmat(BKGR,1,size(SPEC,2));
SPEC(SPEC < 0) = 0;
% Noise estimate and SNR
Noise = mean(mean(abs(SPEC(1:floor(size(SPEC,1)/3),:))));
SNR = 10*log10(SPEC/Noise);
SNR(imag(SNR)~=0) = 0;
SNR(SNR<0) = 0;


%% Strongest Signal

%///////////////////////////////////////////////////////%

SPEC_searchD = SNR;
SPEC_searchD(SPEC_searchD<3) = 0;
[~,IndDop] = max(SPEC_searchD,[],1);

Logo = IndDop > numel(f_lowres)*0.2 & IndDop < numel(f_lowres)*0.8;

for i = 1:length(IndDop)
    if IndDop(i) > length(f_lowres)
        IndDop(i) = length(f_lowres);
    end
end
Sfreq = f_lowres(IndDop);
Measuredfreq = Sfreq(Logo);
time = t(Logo)/86400 + starttimemes;
MeasuredDoppler = Measuredfreq - F_central;
save('saudisatpass4.mat','time','MeasuredDoppler')
