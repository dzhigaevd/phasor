function  [ array xyz] = center_array_com( array,xyonly)
%jclark
%be a little carfeul as the xyz return value is after the centering

if exist('xyonly') ~= 1,xyonly=0;end

array=center_array(array,1); %center according to the max val first

sz=size(array);

nx=sz(1);
dims=ndims(array);

if dims >= 2, ny =sz(2);end
if dims == 3, nz =sz(3);end

dims=ndims(array);

xyz=center_of_mass(abs(array),0);

if mod(nx,2) == 1, xc = (nx+1)/2; else xc =nx/2;end

if dims == 2
    if mod(ny,2) == 1, yc = (ny+1)/2; else yc =ny/2;end
end
if dims == 3
    if mod(ny,2) == 1, yc = (ny+1)/2; else yc =ny/2;end
    if mod(nz,2) ==1, zc = (nz+1)/2; else zc = nz/2;end
end

i=round(xyz(2));
j=round(xyz(1));
if dims == 3,k=round(xyz(3));end
%if dims == 1, array=circshift(array,[-i(1)+xc]);xyz=[-i(1)+xc];end
%if dims == 2, array=circshift(array,[-i(1)+xc,-j(1)+yc]);xyz=[-i(1)+xc,-j(1)+yc];end
%if dims == 3, array=circshift(array,[-i(1)+xc,-j(1)+yc,-k(1)+zc]);xyz=[-i(1)+xc,-j(1)+yc,-k(1)+zc];end
if dims == 1, array=circshift(array,[-i(1)]);xyz=[-i(1)];end
if dims == 2, array=circshift(array,[-i(1),-j(1)]);xyz=[-i(1),-j(1)];end

if xyonly ~= 1
    if dims == 3, array=circshift(array,[-i(1),-j(1),-k(1)]);xyz=[-i(1),-j(1),-k(1)];end
else
    if dims == 3, array=circshift(array,[-i(1),-j(1),0]);xyz=[-i(1),-j(1),-k(1)];end
end
    
%if dims == 3, array=circshift(array,[-i(1),-j(1),-k(1)]);xyz=[-i(1),-j(1),-k(1)];end
   


end

