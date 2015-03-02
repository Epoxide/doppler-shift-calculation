clear all
clc
%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('./functions/dial'));
addpath(genpath('./datasets'));
addpath(genpath('./OrbitCode'));
addpath(genpath('./GPS_CoordinateXforms'));
addpath('./tle');
%%%%%%%%%%%%%%%%%%%%
filename = [cd '/tle/amateur.txt'];
filename2 = [cd '/tle/cubesat.txt'];
urlwrite('http://www.celestrak.com/NORAD/elements/amateur.txt',filename);
urlwrite('http://www.celestrak.com/NORAD/elements/cubesat.txt',filename2);
formatSpec = '%s%[^\n\r]';
% Open the text file.
fileID = fopen(filename,'r');
fileID2 = fopen(filename2,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
dataArray2 = textscan(fileID2, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
% Close the text file.
fclose(fileID);
fclose(fileID2);
% Create output variable
amateur = [dataArray{1:end-1}];
cubesat = [dataArray2{1:end-1}];
n = 1;
for i = 1:3:length(cubesat)
    satNamesCubesat{n,:} = cubesat{i};
    tle1Cubesat{n,:}   	 = cubesat{i+1};
    tle2Cubesat{n,:}   	 = cubesat{i+2};
    n = n + 1;
end
n = 1;
for i = 1:3:length(amateur)
    satNamesAmateur{n,:} = amateur{i};
    tle1Amateur{n,:}     = amateur{i+1};
    tle2Amateur{n,:}     = amateur{i+2};
    n = n + 1;
end
a = [satNamesAmateur; satNamesCubesat];
[satNamesAll,ia,ic] = unique(a,'stable');
satNamesAll = strtrim(satNamesAll);
b = [tle1Amateur;tle1Cubesat];
c = [tle2Amateur;tle2Cubesat];
for i = 1:length(ia)
    tle1All{i,:} = b{ia(i)};
    tle2All{i,:} = c{ia(i)};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TLE Information from Oscar 7

tle1Data = tle1All{1};
tle2Data = tle2All{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
whichconst          = 84;
typerun             = 'c';
typeinput           = 'e';
H                   = 0.500;
mylat               = 59.3496;
mylst               = 18.0724;
Re                  = 6378.137;     % Equatorial Earth's radius [km]
Rp                  = 6356.7523;    % Polar Earth's radius [km]
% Earth's radius at set latitude
Rl                  = sqrt(((Re^2*cosd(mylat))^2+(Rp^2*sind(mylat))^2)/((Re*cosd(mylat))^2+(Rp*sind(mylat))^2))*1e3;
f                   = (Re - Rp)/Re; % Oblateness or flattening
clight              = 299792458;    % Speed of light [m/s]
C1   				= (Re/(1 - (2*f - f^2)*sind(mylat)^2)^0.5 + H)*cosd(mylat);
C2   				= (Re*(1 - f)^2/(1 - (2*f - f^2)*sind(mylat)^2)^0.5 + H)*sind(mylat);
% Position vector of the observer,GEF
R_ob 				= [C1*cosd(mylst), C1*sind(mylst),C2];
% GE_TH is direction cosine matrix to transform position vector components
% from geocentric equatorial frame into the topocentric horizon fream
GE_TH 				= [-sind(mylst)          cosd(mylst)              0;
    -sind(mylat)*cosd(mylst) -sind(mylat)*sind(mylst)  cosd(mylat);
    cosd(mylat)*cosd(mylst)  cosd(mylat)*sind(mylst)   sind(mylat)
    ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% More TLE
[satrec, startmfe, stopmfe, deltamin] = twoline2rv(whichconst, tle1Data,tle2Data,typerun,typeinput);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Azimuth and Elevation
Slrangetable = zeros(1,2);
satname = ['Satellite: ' amateur{1}];
disp(satname)
satfreq = 145.950*1e6; % Satellite frequency [Hz]
for i = 1:2
    time = datetime('now','TimeZone','UTC');
    yr  = time.Year;
    mon = time.Month;
    day = time.Day;
    hr  = time.Hour;
    min = time.Minute;
    sec = time.Second;
    
    jd = 367.0 * yr  ...
        - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 )   ...
        + floor( 275 * mon / 9.0 ) ...
        + day + 1721013.5  ...
        + ( (sec/60.0 + min ) / 60.0 + hr ) / 24.0;
    jdNow = jd;
    tsince = (jdNow-satrec.jdsatepoch)*24*60;
    [satrec, xsat_ecf, vsat_ecf, gst] = spg4_ecf(satrec,tsince);
    R_sc    = xsat_ecf';
    % Position vector of the spacecraft relative to the observer
    R_rel = R_sc - R_ob';
    llhh = ecf2llhT(R_sc'*1e3);
    llh(1) = radtodeg(llhh(1));
    llh(2) = radtodeg(llhh(2));
    R_rel_TH = GE_TH*R_rel;
    rv = R_rel_TH/norm(R_rel_TH);
    Elev = asin(rv(3))*180/pi;      % Elevation angle
    Azf  = atan2(rv(1),rv(2))*180/pi; % Azimuth angle
    Slrange = sqrt((R_rel_TH(1)^2+R_rel_TH(2)^2+R_rel_TH(3)^2))*1e3; % Slant range [m]
    if Azf < 0
        Azf = Azf + 360;
    end
    sat_Alt = -Rl+sqrt(((Rl^2)+(Slrange^2)+(2*Slrange*Rl*sind(Elev)))); % Altitude over Earth [m]
    R_sat   = sqrt(((Re^2*cosd(llh(1)))^2+(Rp^2*sind(llh(1)))^2)/((Re*cosd(llh(1)))^2+(Rp*sind(llh(1)))^2))*1e3;
    sat_Vel = sqrt(3.987*1e14/(R_sat+sat_Alt));
    Slrangetable(i) = Slrange;
    timestring = datestr(time);
    fprintf('\n')
    disp(['Time: ' timestring])
    disp(['Azimuth: ', num2str(Azf) char(176)])
    disp(['Elevation: ', num2str(Elev) char(176)])
    disp(['Slant range: ', num2str(Slrange*1e-3), ' km'])
    disp(['Altitude: ', num2str(sat_Alt*1e-3), ' km'])
    disp(['Velocity: ', num2str(sat_Vel), ' m/s'])
    if i < 2
        pause(10)
    end
end
dSlrange = Slrangetable(1)-Slrangetable(2); % Delta slant range [m]
V_rel = dSlrange/10; % Relative velocity between ground station and satellite [m/s]
fprintf('\n')
disp(['Relative velocity: ' num2str(V_rel), ' m/s']);
Doppler_shift = V_rel*satfreq/clight; % Doppler shift [Hz]
disp(['Predicted Doppler shift: ' num2str(Doppler_shift*1e-3), ' kHz'])
