function [R_dqp_12, R_dqp_3, R_xyz, S_0lab_dir] = plugin_NanoMAX(sampleVerticalAxisRotation, sampleBeamAxisRotation, sampleHorizontalAxisRotation, detectorVerticalAxisRotation, detectorHorizontalAxisRotation, rocking_increment, rocking_angle)
% convert reflection NanoMAX angles to a right-handed coordinate system for new reflection
% adapted by D.Dzhigaev 2020, Lund University
    chi = 90 - sampleBeamAxisRotation;

    % rotation matrices
    R_dqp_12 = rotyd(-detectorVerticalAxisRotation)*rotxd(-detectorHorizontalAxisRotation);
    
    R_xyz = rotzd(chi)*rotyd(-sampleVerticalAxisRotation)*rotxd(-sampleHorizontalAxisRotation); % rotation matrix to rotate a vector in sample coordiantes into lab coordinates
    
    if strcmp(rocking_angle, 'gonphi')
        R_dqp_3 = rotyd(rocking_increment); % it's the negative in rocking increment because we scan from negative to positive
    elseif strcmp(rocking_angle, 'gontheta')
        R_dqp_3 = rotxd(rocking_increment); % it's the negative in rocking increment because we scan from negative to positive
    end

    % X-ray beam incidence direction for NanoMAX: [x,y,z]
    S_0lab_dir = [0; 0; 1];
end