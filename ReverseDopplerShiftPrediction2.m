clear all, clc
%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
amateur2=char(amateur(2));
epochy=amateur2(19:20); %Epoch Year
epochd=amateur2(21:32); %Epoch Day
ftd=amateur2(34:43);
std=amateur2(45:52);
bstar=amateur2(54:61);
amateur3=char(amateur(3));
inclination=amateur3(9:16); %Inclination
rightascension=amateur3(18:25); %Right Ascension of the Ascending Node
eccentricity=amateur3(27:33); %Eccentr




