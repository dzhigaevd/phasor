function meanCrystalSize = scherrerEquation(k,wavelength,FWHM,sampleDetectorDistance,braggAngle)
    meanCrystalSize = k*wavelength/(2*atan(FWHM/(2*sampleDetectorDistance))*cosd(braggAngle));
end

