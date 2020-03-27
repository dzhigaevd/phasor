function [ xyz] = center_of_mass( array,GPU)
%jclark
%array=double(array);
if ~exist('GPU','var')
    GPU = 0;
end

xyz=0;
%MJC edited as part of GPU mod
if GPU==1
 array=gpuArray(array);
 xyz=gpuArray(xyz);
end

tot=sum(sum(sum(array)));

sz=size(array);


nx=sz(2);
ny=sz(1);
if ndims(array) == 3,nz=sz(3);end


if ndims(array) == 3
    [x , y,z]=meshgrid( -(nx-1)/2:(nx-1)/2,-(ny-1)/2:(ny-1)/2,-(nz-1)/2:(nz-1)/2);
end
if ndims(array) == 2
    [x , y]=meshgrid( -(nx-1)/2:(nx-1)/2,-(ny-1)/2:(ny-1)/2);
end

xyz(1)=sum(sum(sum(array.*x)))/tot;
xyz(2)=sum(sum(sum(array.*y)))/tot;

if ndims(array) ==3, xyz(3) = sum(sum(sum(array.*z)))/tot;

end
