function dataTime = getTimeStamp()
%GETTIMESTAMP Summary of this function goes here
%   Detailed explanation goes here
    dataVector = clock;
    dataTime = strcat(num2str(dataVector(3)),'.',num2str(dataVector(2)),'.',num2str(dataVector(1)),'_',num2str(dataVector(4)),'_',num2str(dataVector(5)));    
end