function [ temp ] = pad_by_one(array,drc)
%jclark

switch drc
    case 1
        switch ndims(array)
            case 2
                temp=zeros(size(array)+[1,0]);
                temp(2:end,:)=array(:,:);
            case 3
                temp=zeros(size(array)+[1,0,0]);
                temp(2:end,:,:)=array(:,:,:);
            case 4
                temp=zeros(size(array)+[1,0,0,0]);
                temp(2:end,:,:,:)=array(:,:,:,:);
        end
        
    case 2
        switch ndims(array)
            case 2
                temp=zeros(size(array)+[0,1]);
                temp(:,2:end)=array;
            case 3
                temp=zeros(size(array)+[0,1,0]);
                temp(:,2:end,:)=array;
            case 4
                temp=zeros(size(array)+[0,1,0,0]);
                temp(:,2:end,:,:)=array;
                
        end
        
    case 3
        switch ndims(array)
            case 2
                temp=zeros(size(array)+[0,1]);
                temp(:,2:end)=array;
            case 3
                temp=zeros(size(array)+[0,0,1]);
                temp(:,:,2:end)=array;
            case 4
                temp=zeros(size(array)+[0,0,1,0]);
                temp(:,:,2:end,:)=array;
                
        end
        
    case 4
        switch ndims(array)
            case 2
                temp=zeros(size(array)+[0,1]);
                temp(:,2:end)=array;
            case 3
                temp=zeros(size(array)+[0,0,1]);
                temp(:,:,2:end)=array;
            case 4
                temp=zeros(size(array)+[0,0,0,1]);
                temp(:,:,:,2:end)=array;
                
        end

end

