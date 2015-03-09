function [Doppler_shifts] = ReverseDopplerShiftPrediction(satselect,starttime,endtime,satfreq)
    addpath(genpath('./OrbitCode'));
    addpath(genpath('./GPS_CoordinateXforms'));
    addpath('./tle');
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
    [satNamesAll,ia,~] = unique(a,'stable');
    satNamesAll = strtrim(satNamesAll);
    b = [tle1Amateur;tle1Cubesat];
    c = [tle2Amateur;tle2Cubesat];
    for i = 1:length(ia)
        tle1All{i,:} = b{ia(i)};
        tle2All{i,:} = c{ia(i)};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TLE Information for selected sat
    satnum = strmatch(lower(satselect),lower(satNamesAll));
    while length(satnum) ~= 1
        if length(satnum) > 1
            fprintf('\n')
            disp('Several satellites with that name.')
            satselect = input('Satellite name: ', 's');
            satnum = strmatch(lower(satselect),lower(satNamesAll));
        end
        if isempty(satnum)
            fprintf('\n')
            disp('No satellite with that name.')
            satselect = input('Satellite name: ', 's');
            satnum = strmatch(lower(satselect),lower(satNamesAll));
        end
    end
    fprintf('\n')
    tle1Data = tle1All{satnum};
    tle2Data = tle2All{satnum};

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
    [satrec, ~, ~, ~] = twoline2rv(whichconst, tle1Data,tle2Data,typerun,typeinput);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Azimuth and Elevation

    while length(starttime) ~= 19
        disp('Wrong format')
        starttime = input('Time (YYYY-MM-DD HH:MM:SS in UTC): ','s');
    end
    fprintf('\n')
    while length(endtime) ~= 19
        disp('Wrong format')
        endtime = input('Time (YYYY-MM-DD HH:MM:SS in UTC): ','s');
    end
    fprintf('\n')
    yr  = str2double(starttime(1:4));
    mon = str2double(starttime(6:7));
    day = str2double(starttime(9:10));
    hr  = str2double(starttime(12:13));
    min = str2double(starttime(15:16));
    sec = str2double(starttime(18:19));
    yrend  = str2double(endtime(1:4));
    monend = str2double(endtime(6:7));
    dayend = str2double(endtime(9:10));
    hrend  = str2double(endtime(12:13));
    minend = str2double(endtime(15:16));
    secend = str2double(endtime(18:19));
    if secend < sec && minend < min
        timelength = [hrend-1-hr 60-min-1+minend 60-sec+secend];
        timelength = timelength(1)*3600+timelength(2)*60+timelength(3);
    elseif minend < min
        timelength = [hrend-1-hr 60-min+minend secend-sec];
        timelength = timelength(1)*3600+timelength(2)*60+timelength(3);
    elseif secend < sec
        timelength = [hrend-hr -1+minend-min 60-sec+secend];
        timelength = timelength(1)*3600+timelength(2)*60+timelength(3);
    else
        timelength = [hrend-hr minend-min secend-sec];
        timelength = timelength(1)*3600+timelength(2)*60+timelength(3);
    end
    
    Doppler_shifts = zeros(1,timelength);
    for x = 1:timelength
        jd = 367.0 * yr  ...
            - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 )   ...
            + floor( 275 * mon / 9.0 ) ...
            + day + 1721013.5  ...
            + ( (sec/60.0 + min ) / 60.0 + hr ) / 24.0;
        jdNow = jd;
        tsince = (jdNow-satrec.jdsatepoch)*24*60;
        [~, xsat_ecf, ~, ~] = spg4_ecf(satrec,tsince);
        R_sc    = xsat_ecf';
        % Position vector of the spacecraft relative to the observer
        R_rel = R_sc - R_ob';
        llhh = ecf2llhT(R_sc'*1e3);
        llh(1) = radtodeg(llhh(1));
        llh(2) = radtodeg(llhh(2));
        R_rel_TH = GE_TH*R_rel;
        Slrange = sqrt((R_rel_TH(1)^2+R_rel_TH(2)^2+R_rel_TH(3)^2))*1e3; % Slant range [m]
        
        %%%%%%%%%%%%%%%
        % 1 sec future
        [satrec, ~, ~, ~] = twoline2rv(whichconst, tle1Data,tle2Data,typerun,typeinput);
        if sec == 59
            futuresec = sec-1;
        else
            futuresec = sec+1;
        end
        
        jdfuture = 367.0 * yr  ...
            - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 )   ...
            + floor( 275 * mon / 9.0 ) ...
            + day + 1721013.5  ...
            + ( (futuresec/60.0 + min ) / 60.0 + hr ) / 24.0;
        jdFuture = jdfuture;
        tsinceFuture = (jdFuture-satrec.jdsatepoch)*24*60;
        [~, xsat_ecf, ~, ~] = spg4_ecf(satrec,tsinceFuture);
        R_sc    = xsat_ecf';
        % Position vector of the spacecraft relative to the observer
        R_rel = R_sc - R_ob';
        llhh = ecf2llhT(R_sc'*1e3);
        llh(1) = radtodeg(llhh(1));
        R_rel_TH = GE_TH*R_rel;
        SlrangeFuture = sqrt((R_rel_TH(1)^2+R_rel_TH(2)^2+R_rel_TH(3)^2))*1e3; % Slant range [m]
        
        if sec == 59
            dSlrange = Slrange - SlrangeFuture; % Delta slant range [m]
        else
            dSlrange = SlrangeFuture - Slrange; % Delta slant range [m]
        end
        V_rel = dSlrange; % Relative velocity between ground station and satellite [m/s] (Delta t = 1)
        Doppler_shifts(x) = -V_rel*satfreq*1e6/clight; % Doppler shift [Hz]
        sec = sec + 1;
        if sec == 60
            sec = 00;
            min = min + 1;
            if min == 60
                min = 00;
                hr = hr + 1;
            end
        end
    end
end
