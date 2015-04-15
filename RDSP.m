function [Doppler_shifts] = RDSP(tleData,time,satfreq,OrbitParam)
    addpath(genpath('./OrbitCode'));
    addpath(genpath('./GPS_CoordinateXforms'));
    addpath('./tle');

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
    tle1Data = tleData(1:69);
    tle2Data = tleData(71:end);
    if nargin > 3
        incl = OrbitParam(1);
        if incl >= 100
            tle2Data(9:16) = num2str(incl,'%.4f');
        elseif incl < 100 && incl >= 10
            tle2Data(9) = ' ';
            tle2Data(10:16) = num2str(incl,'%.4f');
        elseif incl < 10
            tle2Data(9:10) = '  ';
            tle2Data(11:16) = num2str(incl,'%.4f');
        end
        
        RightAsc = OrbitParam(2);
        if RightAsc >= 100
            tle2Data(18:25) = num2str(RightAsc,'%.4f');
        elseif RightAsc < 100 && RightAsc >= 10
            tle2Data(18) = ' ';
            tle2Data(19:25) = num2str(RightAsc,'%.4f');
        elseif RightAsc < 10
            tle2Data(18:19) = '  ';
            tle2Data(20:25) = num2str(RightAsc,'%.4f');
        end
        
        ecc = OrbitParam(3);
        ecc = num2str(ecc, '%.7f');
        tle2Data(27:33) = ecc(3:end);
        
        argper = OrbitParam(4);
        if argper >= 100
            tle2Data(35:42) = num2str(argper,'%.4f');
        elseif argper < 100 && argper >= 10
            tle2Data(35) = ' ';
            tle2Data(36:42) = num2str(argper,'%.4f');
        elseif argper < 10
            tle2Data(35:36) = '  ';
            tle2Data(37:42) = num2str(argper,'%.4f');
        end
        
        meanan = OrbitParam(5);
        if meanan >= 100
            tle2Data(44:51) = num2str(meanan,'%.4f');
        elseif meanan < 100 && meanan >= 10
            tle2Data(44) = ' ';
            tle2Data(45:51) = num2str(meanan,'%.4f');
        elseif meanan < 10
            tle2Data(44:45) = '  ';
            tle2Data(46:51) = num2str(meanan,'%.4f');
        end
        
        meanmo = OrbitParam(6);
        if meanmo >= 10
            tle2Data(53:63) = num2str(meanmo,'%.8f');
        elseif meanmo < 10
            tle2Data(53) = ' ';
            tle2Data(54:63) = num2str(meanmo,'%.8f');
        end
    end
    [satrec, ~, ~, ~] = twoline2rv(whichconst, tle1Data,tle2Data,typerun,typeinput);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Doppler_shifts = zeros(1,numel(time));
    for x = 1:numel(time)
        timevect = datevec(time(x));
        yr  = timevect(1);
        mon = timevect(2);
        day = timevect(3);
        hr  = timevect(4);
        min = timevect(5);
        sec = timevect(6);
       
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
        if sec == 59
            futuresec = 00;
            futuremin = min + 1;
            if futuremin == 60
                futuremin = 00;
                futurehr = hr + 1;
            end
        else
            futuresec = sec+1;
            futuremin = min;
            futurehr = hr;
        end
        
        jdfuture = 367.0 * yr  ...
            - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 )   ...
            + floor( 275 * mon / 9.0 ) ...
            + day + 1721013.5  ...
            + ( (futuresec/60.0 + futuremin ) / 60.0 + futurehr ) / 24.0;
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
        dSlrange = SlrangeFuture - Slrange; % Delta slant range [m]
        V_rel = dSlrange; % Relative velocity between ground station and satellite [m/s] (Delta t = 1)
        Doppler_shifts(x) = V_rel*satfreq*1e6/clight; % Doppler shift [Hz]
    end
end
