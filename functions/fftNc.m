	function out = fftNc(in)
		% G = fft2c(g)
		% fft2c performs 2-dimensional Fourier transformation
		% fft2c is normalized (i.e. norm(g) = norm(G) ), 
		% i.e. it preserves the L2-norm
		% if g is two-dimensional, fft2c(g) yields the 2D iDFT of g
		% if g is multi-dimensional, fft2c(g) yields the 2D iDFT of g for each slice
		% along the third dimension
		out = fftshift(fftn(ifftshift(in)));% / sqrt(numel(in));
%         out = fftshift(fftn(ifftshift(in))) / nthroot(numel(in),3);
	end
	