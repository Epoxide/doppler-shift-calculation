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


Sfreq = f_lowres(IndDop);
Measuredfreq = Sfreq(Logo);
time = t(Logo)/86400 + starttimemes;
MeasuredDoppler = Measuredfreq - F_central;
MeasuredDopplerPreFilter = MeasuredDoppler;
for i = 1:length(MeasuredDoppler)
    if MeasuredDoppler(i) > (1e3*12) || MeasuredDoppler(i) < (-1e3*12)
        MeasuredDoppler(i) = NaN;
    end
    if MeasuredDoppler(i) < (200) && MeasuredDoppler(i) > (-200)
        MeasuredDoppler(i) = NaN;
    end
end

Logo = isnan(MeasuredDoppler);
MeasuredDoppler(Logo) = [];
time(Logo) = [];
[tleData, OrbitParam] = Name2TLE('saudisat');
satfreq = 436.795*1e6; % [MHz]
[PredictedDoppler,~,~] = RDSP(tleData, time, satfreq*1e-6, OrbitParam);
MeasuredDopplerOneK = MeasuredDoppler;

for i = 1:length(MeasuredDoppler)
    if abs(MeasuredDoppler(i)-PredictedDoppler(i)) > 1000
        MeasuredDoppler(i) = NaN;
    end
end
Logo = isnan(MeasuredDoppler);
MeasuredDoppler(Logo) = [];
time(Logo) = [];

for i = 2:length(MeasuredDoppler)
    if MeasuredDoppler(i) <= MeasuredDoppler(i-1)
        MeasuredDoppler(i) = NaN;
    end
end
Logo = isnan(MeasuredDoppler);
MeasuredDoppler(Logo) = [];
time(Logo) = [];
[PredictedDoppler, Predsat_Alt, Predsat_Vel] = RDSP(tleData, time, satfreq*1e-6, OrbitParam);
% Fit using lsqcurvefit
beta0 = OrbitParam;
options=optimset('TolFun',1e-8,'TolX',1e-8);
[FittedOrbitParam,RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(@(OrbitParam, time) RDSP(tleData, time, satfreq*1e-6, OrbitParam), beta0, time, MeasuredDoppler, OrbitParam*0.9, OrbitParam*1.1, options);
%FittedOrbitParam = nlinfit(time,MeasuredDoppler,@(OrbitParam, time) RDSP(tleData,time,satfreq*1e-6,OrbitParam), beta0);
[FittedDoppler, Fitsat_Alt, Fitsat_Vel] = RDSP(tleData, time, satfreq*1e-6, FittedOrbitParam);

disp(['Predicted Orbit Parameters:' char(10) 'Inclination: ' num2str(OrbitParam(1)) ' degrees' char(10) 'Right Ascension of the Ascending Node: ' num2str(OrbitParam(2)) ' degrees' char(10) 'Eccentricity: ' num2str(OrbitParam(3)) char(10) 'Argument of Perigee: ' num2str(OrbitParam(4)) ' degrees' char(10) 'Mean Anomaly: ' num2str(OrbitParam(5)) ' degrees' char(10) 'Mean Motion: ' num2str(OrbitParam(6)) ' revs/day']);
disp(' ')
disp(['Fitted Orbit Parameters:' char(10) 'Inclination: ' num2str(FittedOrbitParam(1)) ' degrees' char(10) 'Right Ascension of the Ascending Node: ' num2str(FittedOrbitParam(2)) ' degrees' char(10) 'Eccentricity: ' num2str(FittedOrbitParam(3)) char(10) 'Argument of Perigee: ' num2str(FittedOrbitParam(4)) ' degrees' char(10) 'Mean Anomaly: ' num2str(FittedOrbitParam(5)) ' degrees' char(10) 'Mean Motion: ' num2str(FittedOrbitParam(6)) ' revs/day']);
disp(' ')
Altres = zeros(1,length(Predsat_Alt));
Velres = zeros(1,length(Predsat_Vel));
for i = 1:length(Predsat_Alt)
    Altres(i) = abs(Predsat_Alt(i)-Fitsat_Alt(i));
end
for i = 1:length(Predsat_Vel)
    Velres(i) = abs(Predsat_Vel(i)-Fitsat_Vel(i));
end
Altaverage = mean(Altres);
Velaverage = mean(Velres);
disp(['Average altitude difference from predicted and fitted: ' num2str(Altaverage), ' m']);
disp(['Average velocity difference from predicted and fitted: ' num2str(Velaverage), ' km/s']);

figure(1)
hold on
plot(time,MeasuredDoppler,'rx',time,PredictedDoppler,'bx',time,FittedDoppler,'gx')
legend('Measured Doppler Shift','Predicted Doppler Shift','Fitted Doppler Shift')
datetick('x', 13)
title('Measured, Predicted, and Fitted Doppler Shifts')
xlabel('Time')
ylabel('Doppler Shift [Hz]')
figure(2)
hold on
subplot(2,1,1)
plot(time,Predsat_Alt,'rx',time,Fitsat_Alt,'bx')
legend('Predicted Satellite Altitude','Fitted Satellite Altitude')
datetick('x', 13)
title('Predicted, and Fitted Satellite Altitudes')
xlabel('Time')
ylabel('Satellite Altitude [m]')
subplot(2,1,2)
plot(time,Predsat_Vel,'rx',time,Fitsat_Vel,'bx')
legend('Predicted Satellite Orbital Velocity','Fitted Satellite Orbital Velocity Velocity')
datetick('x', 13)
title('Predicted, and Fitted Satellite Orbital Velocities')
xlabel('Time')
ylabel('Satellite Orbital Velocity [km/s]')
