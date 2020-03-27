function savemat2vtk(filename,array,spacing,array2)
%  adapted by DDzhigaev for PHASOR.
%  now accept additional matrix to save (e.g. displacement, strain)
%  jclark
%  adapted from savevtk.m, no does two scalars and set spacing
%  savevtk Save a 3-D scalar array in VTK format.
%  savevtk(array, filename) saves a 3-D array of any size to
%  filename in VTK format.

    if exist('spacing','var') ~= 1
        spacing=1;
    end
    
    disp('Writing Amplitude to VTK file....')
    disp(['Spacing - ',num2str(spacing)])

    [ny, nx, nz] = size(array);
   
    fid = fopen(filename, 'wt');
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Comment goes here\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, '\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
    fprintf(fid, '\n');
    fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
    %fprintf(fid, 'SPACING    1.000   1.000   1.000\n');
    fprintf(fid, 'SPACING    %d   %d   %d\n',spacing,spacing,spacing);
    fprintf(fid, '\n');
    fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
    fprintf(fid, 'SCALARS amp double\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '\n');
    
    for a=1:nz        
            fprintf(fid, '%d ', abs(permute(array(:,:,a),[2 1])));
            %fwrite(fid,permute(array(:,:,a),[2 1]),'double');
            fprintf(fid, '\n');
    end
    
    if ~isreal(array)                   
            disp('Writing Phase to VTK file....')            
            fprintf(fid, 'FIELD FieldData 1 \n');
            fprintf(fid, 'phases 1 %d', nx*ny*nz);
            fprintf(fid, 'double \n');
            for a=1:nz
                    fprintf(fid, '%d ', angle(permute(array(:,:,a),[2 1])));
                    %fwrite(fid, permute(array2(:,:,a),[2 1]),'double');
                    fprintf(fid, '\n');
            end            
    end    
    
    if exist('array2','var') 
            disp('Writing additional value to VTK file....')            
            fprintf(fid, 'FIELD FieldData 1 \n');
            fprintf(fid, 'values 1 %d', nx*ny*nz);
            fprintf(fid, 'double \n');
            for a=1:nz
                    fprintf(fid, '%d ', permute(array2(:,:,a),[2 1]));
                    %fwrite(fid, permute(array2(:,:,a),[2 1]),'double');
                    fprintf(fid, '\n');
            end
    end
    fclose(fid);
    
return

