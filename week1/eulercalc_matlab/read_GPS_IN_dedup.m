function [lon,lat,uE,uN,sE,sN,corr,name]=read_GPS_IN_dedup(filename)

delimiter = ' ';
formatSpec = '%f%f%f%f%f%f%f%s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

lon = dataArray{:, 1};
lat = dataArray{:, 2};
uE = dataArray{:, 3};
uN = dataArray{:, 4};
sE = dataArray{:, 5};
sN = dataArray{:, 6};
corr = dataArray{:, 7};
name = dataArray{:, 8};


%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;
end