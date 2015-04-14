function  [tleData,OrbitParam] = Name2TLE(satselect)
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
    tleData = [tle1Data ' ' tle2Data];
    incl = str2double(tle2Data(9:16));
    RightAsc = str2double(tle2Data(18:25));
    ecc = str2double(tle2Data(27:33))*10^-7;
    argper = str2double(tle2Data(35:42));
    meanan = str2double(tle2Data(44:51));
    meanmo = str2double(tle2Data(53:63));
    OrbitParam = [incl RightAsc ecc argper meanan meanmo];
end

