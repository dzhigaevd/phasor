
    function out = ifft2c(in)
        if mod(size(in,1),2)==0 && mod(size(in,2),2)==0
            % for an even number of entries the 2d inverse FT can be calculated by
            % means of the 2d forward FT (slightly faster)
            out = conj(fftshift(fft2(ifftshift(conj(in))))) / sqrt(numel(in));
        else
            out = fftshift(ifft2(ifftshift(in))) * sqrt(numel(in));
        end
    end  