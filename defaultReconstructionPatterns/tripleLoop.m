%% Simple Phase retrieval pattern
tic
sigma = linspace(5., 3., 5 );

for ii = 1:numel(sigma)     
    ph.ER( 30, 1 )
    ph.SW( sigma(ii), 0.1 )     
end

ph.HIO(300,0.9,1)  

sigma = linspace( 3., 2., 4 );
for ii = 1:numel(sigma)                  
    ph.SF( 25, 1 );
    ph.SW( sigma(ii), 0.1 );      
end

ph.HIO( 300, 0.9, 1)  

sigma = linspace( 2., 1., 5 );
for ii = 1:numel(sigma)                  
    ph.ER( 90, 1 )
    ph.SW( sigma(ii), 0.1 )  
end
toc
disp('Phasing is finished!')