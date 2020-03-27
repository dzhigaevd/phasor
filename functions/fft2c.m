% Complex-valued normalized Fourier transform functions
	function out = fft2c(in)
		% fft2c performs 2-dimensional Fourier transformation
		% fft2c is normalized (i.e. norm(g) = norm(G) ), 
		% i.e. it preserves the L2-norm
		% if g is two-dimensional, fft2c(g) yields the 2D iDFT of g
		% if g is multi-dimensional, fft2c(g) yields the 2D iDFT of g for each slice
		% along the third dimension
		out = fftshift(fft2(ifftshift(in))) / sqrt(numel(in));
	end