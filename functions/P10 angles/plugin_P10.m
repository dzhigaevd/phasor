function [R_dqp_12, R_dqp_3, R_xyz, S_0lab_dir] = plugin_P10(sampleVerticalAxisRotation, sampleBeamAxisRotation, sampleHorizontalAxisRotation, detectorVerticalAxisRotation, detectorHorizontalAxisRotation, rocking_increment, rocking_angle)
    % convert reflection SPEC angles to a right-handed coordinate system for new reflection
    chi = 90-sampleBeamAxisRotation;

    % rotation matrices
    R_dqp_12 = rotyd(detectorVerticalAxisRotation)*rotxd(-detectorHorizontalAxisRotation);
    
    R_xyz = rotyd(sampleVerticalAxisRotation)*rotzd(chi)*rotxd(-sampleHorizontalAxisRotation); % rotation matrix to rotate a vector in sample coordiantes into lab coordinates
    if strcmp(rocking_angle, 'phi')
        R_dqp_3 = rotyd(-rocking_increment); % it's the negative in rocking increment because we scan from negative to positive
    elseif strcmp(rocking_angle, 'om')
        R_dqp_3 = rotxd(rocking_increment); % it's the negative in rocking increment because we scan from negative to positive
    end

    % S_0lab direction for APS
    S_0lab_dir = [0; 0; 1];
end