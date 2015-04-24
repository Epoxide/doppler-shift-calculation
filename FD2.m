clear all, close all, clc
load saudisatpass1.mat
time1 = time;
MeasuredDoppler1 = MeasuredDoppler;
load saudisatpass2.mat
time2 = time;
MeasuredDoppler2 = MeasuredDoppler;
load saudisatpass3.mat
time3 = time;
MeasuredDoppler3 = MeasuredDoppler;
load saudisatpass4.mat
time4 = time;
MeasuredDoppler4 = MeasuredDoppler;

satfreq = 436.795*1e6; % [MHz]

for i = 1:length(MeasuredDoppler1)
    if MeasuredDoppler1(i) > (1e3*12) || MeasuredDoppler1(i) < (-1e3*12)
        MeasuredDoppler1(i) = NaN;
    end
    if MeasuredDoppler1(i) < (200) && MeasuredDoppler1(i) > (-200)
        MeasuredDoppler1(i) = NaN;
    end
end

Logo = isnan(MeasuredDoppler1);
MeasuredDoppler1(Logo) = [];
time1(Logo) = [];

for i = 1:length(MeasuredDoppler2)
    if MeasuredDoppler2(i) > (1e3*12) || MeasuredDoppler2(i) < (-1e3*12)
        MeasuredDoppler2(i) = NaN;
    end
    if MeasuredDoppler2(i) < (200) && MeasuredDoppler2(i) > (-200)
        MeasuredDoppler2(i) = NaN;
    end
end

Logo = isnan(MeasuredDoppler2);
MeasuredDoppler2(Logo) = [];
time2(Logo) = [];

for i = 1:length(MeasuredDoppler3)
    if MeasuredDoppler3(i) > (1e3*12) || MeasuredDoppler3(i) < (-1e3*12)
        MeasuredDoppler3(i) = NaN;
    end
    if MeasuredDoppler3(i) < (200) && MeasuredDoppler3(i) > (-200)
        MeasuredDoppler3(i) = NaN;
    end
end

Logo = isnan(MeasuredDoppler3);
MeasuredDoppler3(Logo) = [];
time3(Logo) = [];

for i = 1:length(MeasuredDoppler4)
    if MeasuredDoppler4(i) > (1e3*12) || MeasuredDoppler4(i) < (-1e3*12)
        MeasuredDoppler4(i) = NaN;
    end
    if MeasuredDoppler4(i) < (200) && MeasuredDoppler4(i) > (-200)
        MeasuredDoppler4(i) = NaN;
    end
end

Logo = isnan(MeasuredDoppler4);
MeasuredDoppler4(Logo) = [];
time4(Logo) = [];

[tleData, OrbitParam] = Name2TLE('saudisat');

[PredictedDoppler1,~,~] = RDSP(tleData, time1, satfreq*1e-6, OrbitParam);
[PredictedDoppler2,~,~] = RDSP(tleData, time2, satfreq*1e-6, OrbitParam);
[PredictedDoppler3,~,~] = RDSP(tleData, time3, satfreq*1e-6, OrbitParam);
[PredictedDoppler4,~,~] = RDSP(tleData, time4, satfreq*1e-6, OrbitParam);

for i = 1:length(MeasuredDoppler1)
    if abs(MeasuredDoppler1(i)-PredictedDoppler1(i)) > 1000
        MeasuredDoppler1(i) = NaN;
    end
end
Logo = isnan(MeasuredDoppler1);
MeasuredDoppler1(Logo) = [];
time1(Logo) = [];

for i = 2:length(MeasuredDoppler1)
    if MeasuredDoppler1(i) < MeasuredDoppler1(i-1)
        MeasuredDoppler1(i) = NaN;
    end
end
Logo = isnan(MeasuredDoppler1);
MeasuredDoppler1(Logo) = [];
time1(Logo) = [];
[PredictedDoppler1, Predsat_Alt1, Predsat_Vel1] = RDSP(tleData, time1, satfreq*1e-6, OrbitParam);

for i = 1:length(MeasuredDoppler2)
    if abs(MeasuredDoppler2(i)-PredictedDoppler2(i)) > 1000
        MeasuredDoppler2(i) = NaN;
    end
end
Logo = isnan(MeasuredDoppler2);
MeasuredDoppler2(Logo) = [];
time2(Logo) = [];

for i = 2:length(MeasuredDoppler2)
    if MeasuredDoppler2(i) < MeasuredDoppler2(i-1)
        MeasuredDoppler2(i) = NaN;
    end
end
Logo = isnan(MeasuredDoppler2);
MeasuredDoppler2(Logo) = [];
time2(Logo) = [];
[PredictedDoppler2, Predsat_Alt2, Predsat_Vel2] = RDSP(tleData, time2, satfreq*1e-6, OrbitParam);

for i = 1:length(MeasuredDoppler3)
    if abs(MeasuredDoppler3(i)-PredictedDoppler3(i)) > 1000
        MeasuredDoppler3(i) = NaN;
    end
end
Logo = isnan(MeasuredDoppler3);
MeasuredDoppler3(Logo) = [];
time3(Logo) = [];

for i = 2:length(MeasuredDoppler3)
    if MeasuredDoppler3(i) < MeasuredDoppler3(i-1)
        MeasuredDoppler3(i) = NaN;
    end
end
Logo = isnan(MeasuredDoppler3);
MeasuredDoppler3(Logo) = [];
time3(Logo) = [];
[PredictedDoppler3, Predsat_Alt3, Predsat_Vel3] = RDSP(tleData, time3, satfreq*1e-6, OrbitParam);

for i = 1:length(MeasuredDoppler4)
    if abs(MeasuredDoppler4(i)-PredictedDoppler4(i)) > 1000
        MeasuredDoppler4(i) = NaN;
    end
end
Logo = isnan(MeasuredDoppler4);
MeasuredDoppler4(Logo) = [];
time4(Logo) = [];

for i = 2:length(MeasuredDoppler4)
    if MeasuredDoppler4(i) < MeasuredDoppler4(i-1)
        MeasuredDoppler4(i) = NaN;
    end
end
Logo = isnan(MeasuredDoppler4);
MeasuredDoppler4(Logo) = [];
time4(Logo) = [];
[PredictedDoppler4, Predsat_Alt4, Predsat_Vel4] = RDSP(tleData, time4, satfreq*1e-6, OrbitParam);

time = [time1, time2, time3, time4];
MeasuredDoppler = [MeasuredDoppler1, MeasuredDoppler2, MeasuredDoppler3, MeasuredDoppler4];
PredictedDoppler = [PredictedDoppler1, PredictedDoppler2, PredictedDoppler3, PredictedDoppler4];

% Fit using lsqcurvefit
beta0 = [1+(359-1).*rand(1) 1+(359-1).*rand(1) 0.000001+(0.1-0.000001).*rand(1) 1+(359-1).*rand(1) 1+(359-1).*rand(1) 9+(20-9).*rand(1)];
for i = 1:2
    while beta0(i) > OrbitParam(i)*0.9 && beta0(i) < OrbitParam(i)*1.1
        beta0(i) = 1+(359-1).*rand(1);
    end
end
while beta0(3) > OrbitParam(3)*0.9 && beta0(3) < OrbitParam(3)*1.1
    beta0(3) = 0.000001+(0.1-0.000001).*rand(1);
end
for i = 4:5
    while beta0(i) > OrbitParam(i)*0.9 && beta0(i) < OrbitParam(i)*1.1
        beta0(i) = 1+(359-1).*rand(1);
    end
end
while beta0(6) > OrbitParam(6)*0.9 && beta0(6) < OrbitParam(6)*1.1
    beta0(6) = 9+(20-9).*rand(1);
end
options=optimset('TolFun',1e-12,'TolX',1e-12);
% lsqcurvefit
[FittedOrbitParam,RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(@(OrbitParam, time) RDSP(tleData, time, satfreq*1e-6, OrbitParam), beta0, time, MeasuredDoppler, OrbitParam*0.9, OrbitParam*1.1, options);
[FittedDoppler, Fitsat_Alt, Fitsat_Vel] = RDSP(tleData, time, satfreq*1e-6, FittedOrbitParam);
[betaDoppler, betasat_Alt, betasat_Vel] = RDSP(tleData, time, satfreq*1e-6, beta0);

disp(['Predicted Orbit Parameters:' char(10) 'Inclination: ' num2str(beta0(1)) ' degrees' char(10) 'Right Ascension of the Ascending Node: ' num2str(beta0(2)) ' degrees' char(10) 'Eccentricity: ' num2str(beta0(3)) char(10) 'Argument of Perigee: ' num2str(beta0(4)) ' degrees' char(10) 'Mean Anomaly: ' num2str(beta0(5)) ' degrees' char(10) 'Mean Motion: ' num2str(beta0(6)) ' revs/day']);
disp(' ');
disp(['Two-line Element Orbit Parameters:' char(10) 'Inclination: ' num2str(OrbitParam(1)) ' degrees' char(10) 'Right Ascension of the Ascending Node: ' num2str(OrbitParam(2)) ' degrees' char(10) 'Eccentricity: ' num2str(OrbitParam(3)) char(10) 'Argument of Perigee: ' num2str(OrbitParam(4)) ' degrees' char(10) 'Mean Anomaly: ' num2str(OrbitParam(5)) ' degrees' char(10) 'Mean Motion: ' num2str(OrbitParam(6)) ' revs/day']);
disp(' ');
disp(['Fitted Orbit Parameters:' char(10) 'Inclination: ' num2str(FittedOrbitParam(1)) ' degrees' char(10) 'Right Ascension of the Ascending Node: ' num2str(FittedOrbitParam(2)) ' degrees' char(10) 'Eccentricity: ' num2str(FittedOrbitParam(3)) char(10) 'Argument of Perigee: ' num2str(FittedOrbitParam(4)) ' degrees' char(10) 'Mean Anomaly: ' num2str(FittedOrbitParam(5)) ' degrees' char(10) 'Mean Motion: ' num2str(FittedOrbitParam(6)) ' revs/day']);

subplot(2,2,1)
plot(time,MeasuredDoppler,'rx',time,PredictedDoppler,'bx',time,FittedDoppler,'gx',time,betaDoppler,'x')
legend('Measured Doppler Shift','Two-line Element Doppler Shift','Fitted Doppler Shift','Initial Guess', 'Location', 'NorthWest')
date1 = datestr(time1(1));
date1 = date1(1:11);
title1 = ['First pass on ' date1];
title(title1)
xlabel('Time')
ylabel('Doppler Shift [Hz]')
xlim([time1(1) time1(end)])
datetick('x', 13, 'keeplimits','keepticks')
subplot(2,2,2)
plot(time,MeasuredDoppler,'rx',time,PredictedDoppler,'bx',time,FittedDoppler,'gx',time,betaDoppler,'x')
legend('Measured Doppler Shift','Two-line Element Doppler Shift','Fitted Doppler Shift','Initial Guess', 'Location', 'NorthWest')
date2 = datestr(time2(1));
date2 = date2(1:11);
title2 = ['Second pass on ' date2];
title(title2)
xlabel('Time')
ylabel('Doppler Shift [Hz]')
xlim([time2(1) time2(end)])
datetick('x', 13, 'keeplimits','keepticks')
subplot(2,2,3)
plot(time,MeasuredDoppler,'rx',time,PredictedDoppler,'bx',time,FittedDoppler,'gx',time,betaDoppler,'x')
legend('Measured Doppler Shift','Two-line Element Doppler Shift','Fitted Doppler Shift','Initial Guess', 'Location', 'NorthWest')
date3 = datestr(time3(1));
date3 = date3(1:11);
title3 = ['Third pass on ' date3];
title(title3)
xlabel('Time')
ylabel('Doppler Shift [Hz]')
xlim([time3(1) time3(end)])
datetick('x', 13, 'keeplimits','keepticks')
subplot(2,2,4)
plot(time,MeasuredDoppler,'rx',time,PredictedDoppler,'bx',time,FittedDoppler,'gx',time,betaDoppler,'x')
legend('Measured Doppler Shift','Two-line Element Doppler Shift','Fitted Doppler Shift','Initial Guess', 'Location', 'NorthWest')
date4 = datestr(time4(1));
date4 = date4(1:11);
title4 = ['Fourth pass on ' date4];
title(title4)
suptitle('Measured, Two-line Element, and Fitted Doppler Shifts with Initial Guess for 4 passes of Saudisat 1C')
xlabel('Time')
ylabel('Doppler Shift [Hz]')
xlim([time4(1) time4(end)])
datetick('x', 13, 'keeplimits','keepticks')
