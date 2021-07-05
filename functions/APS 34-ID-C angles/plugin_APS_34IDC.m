function [R_dqp_12, R_dqp_3, R_xyz, S_0lab_dir] = plugin_APS_34IDC(sampleVerticalAxisRotation, sampleBeamAxisRotation, sampleHorizontalAxisRotation, detectorVerticalAxisRotation, detectorHorizontalAxisRotation, rocking_increment, rocking_angle)
    
    % convert reflection SPEC angles to a right-handed coordinate system for new reflection
    chi = 90 - sampleBeamAxisRotation;

%     % rotation matrices
    R_dqp_12 = rotyd(detectorVerticalAxisRotation)*rotxd(-detectorHorizontalAxisRotation);
    
    R_xyz = rotzd(chi)*rotyd(sampleHorizontalAxisRotation)*rotxd(sampleVerticalAxisRotation); % rotation matrix to rotate a vector in sample coordiantes into lab coordinates
    if strcmp(rocking_angle, 'dtheta')
        R_dqp_3 = rotyd(-rocking_increment); % it's the negative in rocking increment because we scan from negative to positive
    elseif strcmp(rocking_angle, 'dphi')
        R_dqp_3 = rotxd(-rocking_increment); % it's the negative in rocking increment because we scan from negative to positive
    end

    % S_0lab direction for APS
    S_0lab_dir = [0; 0; 1];

end