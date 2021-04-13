function out = readMultiTiff(fname)
    info = imfinfo(fname);
    num_images = numel(info);
    for k = 1:num_images
        out(:,:,k) = imread(fname, k, 'Info', info);        
    end
end