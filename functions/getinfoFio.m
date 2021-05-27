function [ fioParam ] = getinfoFio( fioPath )
%   GETINFO_FIO Summary of this function goes here
%   Detailed explanation goes here
    [data, counters, cstr, fioParam.motorNames, fioParam.motorVals] = readfio(fioPath);
    
%     fioParam.motorPos = scanParam;
    
    fastMotor_indf  = 1;
    slowMotor_indf  = 2;

    try
        fluoDetector_indf = find(ismember(counters,fluoDetector) == 1);
    catch
        disp('Warning! Desired fluorescence detector is not found!')
    end    

    fioParam.fastMotorName   = char(counters{fastMotor_indf});
    fioParam.slowMotorName   = char(counters{slowMotor_indf});

    commandList     = strsplit(cstr);
    fioParam.scanType    = commandList{1};
    disp(['Processing: ',cstr]);

    switch fioParam.scanType
        case {'amesh', 'dmesh'}   
            fioParam.scanType1       = 'mesh';
            fioParam.fastMotor       = commandList{2};
            fioParam.fastMotorStart  = commandList{3};
            fioParam.fastMotorEnd    = commandList{4};
            fioParam.fastMotorPoints = str2num(commandList{5})+1;

            fioParam.slowMotor       = commandList{6};
            fioParam.slowMotorStart  = commandList{7};
            fioParam.slowMotorEnd    = commandList{8};                
            fioParam.slowMotorPoints = str2num(commandList{9})+1;
            
            fioParam.nPatterns       = fioParam.slowMotorPoints*fioParam.fastMotorPoints;
            
            fioParam.coordinates = reshape(data(:,1:2),[fioParam.fastMotorPoints, fioParam.slowMotorPoints, 2]);

            fioParam.xVector = fioParam.coordinates(:,1,1)';
            fioParam.yVector = fioParam.coordinates(1,:,2);
            
            if strcmp(fioParam.fastMotor(end),'y')
                fioParam.nH = fioParam.fastMotorPoints; % horizontal points for linescan checks
                fioParam.nV = fioParam.slowMotorPoints; % vertical points for linescan checks
            else
                fioParam.nH = fioParam.slowMotorPoints; % horizontal points for linescan checks
                fioParam.nV = fioParam.fastMotorPoints; % vertical points for linescan checks
            end                       
        
        case {'ascan', 'dscan', 'a2scan', 'd2scan'}            
            fioParam.scanType1    = 'scan';
                                      
            switch fioParam.scanType 
                case {'ascan', 'dscan'} 
                    fioParam.Motor       = commandList{2};
                    fioParam.MotorStart  = commandList{3};
                    fioParam.MotorEnd    = commandList{4};
                    fioParam.nPatterns   = str2num(commandList{5})+1; % vertical points for linescan checks  
                    fioParam.Vector      = data(:,1);
                    
                case {'a2scan', 'd2scan'}
                    fioParam.nPatterns   = str2num(commandList{8})+1; % vertical points for linescan checks 
                    fioParam.Vector1      = data(:,1);
                    fioParam.Vector2      = data(:,2);
            end
      
        case {'series'}
            fioParam.nP = str2num(commandList{5})+1; % vertical points for linescan checks   
            
        otherwise
            disp('Unknown scantype');
    end    
end

