


% params.lam = postPhasor.experiment.wavelength;
% params.dtilt = 0;
% params.det_px = 55*10^-6;
% params.gam = params.gamma;
% params.binning = [1 1 1];
% params.GPU = 1;
% params.arm = params.arm/1000;
%[pn params] = det2lab(pn,params);

Nx = size(postPhasor.object,1);
Ny = size(postPhasor.object,2);
Nz = size(postPhasor.object,3);

% lambda = params.lam;         %wavelength (m)
% arm = params.arm;                      %camera distance (m)
% det_pix = params.det_px;               %detector pixel size (m)
% gamma = params.gam;                     %gamma (degrees)
% delta = params.delta;                  %delta (degrees)
% dth = params.dth;                      %angular step size
dtilt = 0;

px = postPhasor.plotting.binning(1) * postPhasor.experiment.detector_pitch;
py = postPhasor.plotting.binning(2) * postPhasor.experiment.detector_pitch;

sx=postPhasor.experiment.wavelength*postPhasor.experiment.sample_detector_d/px/Nx;
sy=postPhasor.experiment.wavelength*postPhasor.experiment.sample_detector_d/py/Ny;

if postPhasor.experiment.angular_step ~= 0,
    sz=abs(postPhasor.experiment.wavelength/(postPhasor.experiment.angular_step*pi/180)/Nz);
else
    sz=abs(postPhasor.experiment.wavelength/(dtilt*pi/180)/Nz);
end

FOVth=0.66;      %maximum reduction allowed in the field of view to prevent cutting the object
sample_pixel=min([sz,sy,sx]);
FOVold=[Nx*sx,Ny*sy,Nz*sz];
FOVnew=[Nx*sample_pixel,Ny*sample_pixel,Nz*sample_pixel];
FOVratio=min(FOVnew./FOVold);

if FOVratio(1) < FOVth,
    sample_pixel=FOVth/FOVratio(1)*sample_pixel;
end

dx=1/Nx;
dy=1/Ny;
dz=1/Nz;

deg2rad=pi/180;

postPhasor.experiment.delta=postPhasor.experiment.delta*deg2rad;
postPhasor.experiment.gamma=postPhasor.experiment.gamma*deg2rad;
postPhasor.experiment.angular_step=postPhasor.experiment.angular_step*deg2rad;

dpx=px/postPhasor.experiment.sample_detector_d;
dpy=py/postPhasor.experiment.sample_detector_d;

%old pre july 2013
dQdpx(1) = -cos(postPhasor.experiment.delta)*cos(postPhasor.experiment.gamma);
dQdpx(2) = 0.0;
dQdpx(3) = +sin(postPhasor.experiment.delta)*cos(postPhasor.experiment.gamma);

dQdpy(1) = sin(postPhasor.experiment.delta)*sin(postPhasor.experiment.gamma);
dQdpy(2) = -cos(postPhasor.experiment.gamma);
dQdpy(3) = cos(postPhasor.experiment.delta)*sin(postPhasor.experiment.gamma);

dQdth(1) = -cos(postPhasor.experiment.delta)*cos(postPhasor.experiment.gamma)+1.0;
dQdth(2) = 0.0;
dQdth(3) = sin(postPhasor.experiment.delta)*cos(postPhasor.experiment.gamma);
 
Astar(1) = (2*pi/postPhasor.experiment.wavelength)*dpx*dQdpx(1);
Astar(2) = (2*pi/postPhasor.experiment.wavelength)*dpx*dQdpx(2);
Astar(3) = (2*pi/postPhasor.experiment.wavelength)*dpx*dQdpx(3);

Bstar(1) = (2*pi/postPhasor.experiment.wavelength)*dpy*dQdpy(1);
Bstar(2) = (2*pi/postPhasor.experiment.wavelength)*dpy*dQdpy(2);
Bstar(3) = (2*pi/postPhasor.experiment.wavelength)*dpy*dQdpy(3);

Cstar(1) = (2*pi/postPhasor.experiment.wavelength)*postPhasor.experiment.angular_step*dQdth(1);
Cstar(2) = (2*pi/postPhasor.experiment.wavelength)*postPhasor.experiment.angular_step*dQdth(2);
Cstar(3) = (2*pi/postPhasor.experiment.wavelength)*postPhasor.experiment.angular_step*dQdth(3);

denom=dot(Astar,cross(Bstar,Cstar));
Axdenom=cross(Bstar,Cstar);
Bxdenom=cross(Cstar,Astar);
Cxdenom=cross(Astar,Bstar);

A(1)=2*pi*Axdenom(1)/(denom);
A(2)=2*pi*Axdenom(2)/(denom);
A(3)=2*pi*Axdenom(3)/(denom);

B(1)=2*pi*Bxdenom(1)/(denom);
B(2)=2*pi*Bxdenom(2)/(denom);
B(3)=2*pi*Bxdenom(3)/(denom);

C(1)=2*pi*Cxdenom(1)/(denom);
C(2)=2*pi*Cxdenom(2)/(denom);
C(3)=2*pi*Cxdenom(3)/(denom);

T1=[A(1) B(1) C(1) 0;A(2) B(2) C(2) 0;A(3) B(3) C(3) 0;0 0 0 1];

T = T1(1:3,1:3);

dx = 1/Nx;
dy = 1/Ny;
dz = 1/Nz;


[x y z]=meshgrid( ((1:Nx)*dx),((1:Ny)*dy),(1:Nz)*dz);

r = zeros([3, Nx, Ny, Nz]);

r(1,:,:,:) = x;
r(2,:,:,:) = y;
r(3,:,:,:) = z;

r = reshape(r, 3, Nx*Ny*Nz);

coords = T*r;

coords = reshape(coords, 3, Nx, Ny, Nz);

X = reshape(coords(1,:), Nx, Ny, Nz);
Y = reshape(coords(2,:), Nx, Ny, Nz);
Z = reshape(coords(3,:), Nx, Ny, Nz);

X = X - min(X(:));
Y = Y - min(Y(:));
Z = Z - min(Z(:));

XX=x*sample_pixel*Nx;
YY=y*sample_pixel*Ny;
ZZ=z*sample_pixel*Nz;

postPhasor.object=flipdim(postPhasor.object,3);

v1=X-max(X(:))/2;
v2=Y-max(Y(:))/2;
v3=Z-max(Z(:))/2;
v4=XX-max(XX(:))/2;
v5=YY-max(YY(:))/2;
v6=ZZ-max(ZZ(:))/2;
%%

%%
w = griddata(v1,v2,v3,postPhasor.object,v4,v5,v6);

w(isnan(w))=0;

%%
isoval = 0.2;
amp = abs(w);
    ph = angle(w);  
    sz = size(w);
    %[X Y Z]=meshgrid(pixsize*[-sz(1)/2:sz(1)/2-1], pixsize*[-sz(2)/2:sz(2)/2-1], pixsize*[-sz(3)/2:sz(3)/2-1]);
    s1 = subplot(1,1,1)
    h1 = patch(isosurface(amp/max(amp(:)), isoval));
    isonormals(amp, h1);
    isocolors(ph, h1);
    set(h1, 'facecolor', 'interp');
    set(h1, 'edgecolor', 'none');
    axis image;
    camlight right;
    %caxis([-pi pi]);
    lighting gouraud
    set(h1, 'SpecularStrength', .1);
    camlight
   % title(['Phases (', num2str(isoval), '), Scan N', num2str(filenums(ii))]);
    set(gca,'Fontsize',18)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    colorbar
    colormap(parula)
    %myColorMap = jet(256);
    %myColorMap(1,:) = [0.8 0.8 0.8];
    %colormap(s1, myColorMap);
    
    %Create light
light('Position',[1797.82460021973 0 -692.443905290251],'Style','local');

%Create light
light('Position',[-1795.06953430176 0 -692.443905290252],'Style','local');

%Create light
light('Position',[1.37753295898438 2598.07621135332 -6.08567047119132],...
   'Style','local');

%Create light  
light('Position',[1.37753295898438 -2598.07621135332 -6.08567047119132],...
   'Style','local');
 view([-27, -17])   
 
  view([0, 89])   