   function out = ifftNc(in)
        if mod(size(in,1),2)==0 && mod(size(in,2),2)==0 && mod(size(in,3),2)==0
            % for an even number of entries the 2d inverse FT can be calculated by
            % means of the 2d forward FT (slightly faster)
            out = conj(fftshift(fftn(ifftshift(conj(in))))) / sqrt(numel(in));
%             out = conj(fftshift(fftn(ifftshift(conj(in))))) / nthroot(numel(in),3);
        else
            out = fftshift(ifftn(ifftshift(in)))* sqrt(numel(in));
%             out = fftshift(ifftn(ifftshift(in))) * nthroot(numel(in),3);
        end
    end