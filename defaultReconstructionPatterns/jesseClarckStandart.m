%% Simple Phase retrieval pattern
tic

iteration = 620;
threshold = 0.1;
nER = 20;
nHIO = 180;
nIteration = round(iteration/(2*nER+nHIO))

sigma = linspace(5., 1., nIteration*4 );
sigma_inc = 1;

for ii = 1:nIteration 
    fprintf('Overall iteration: %d\n',ii)
    ph.ER(nER, 0)
    ph.SW( sigma(sigma_inc), threshold )
    sigma_inc = sigma_inc+1;
    
    ph.HIO(nHIO, 0)
    ph.SW( sigma(sigma_inc), threshold/2 ) 
    sigma_inc = sigma_inc+1;
    
    ph.SF(nER, 0)
    ph.SW( sigma(sigma_inc), threshold )
    sigma_inc = sigma_inc+1;
    
    ph.HIO(nHIO, 0)
    ph.SW( sigma(sigma_inc), threshold/2 ) 
    sigma_inc = sigma_inc+1;
end

toc
disp('Phasing is finished!')