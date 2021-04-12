function output = crop3D(input)
% for convex objects
    input_amp = abs(input);

    figure; 
    imagesc(squeeze(sum(input_amp,3)));
    
    h = imrect;
    
    pos = round(getPosition(h)); %[xmin ymin width height]
    
    for ii = 1:size(input,3)
        output(:,:,ii) = squeeze(input(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3),ii));
    end
    close;        
end

