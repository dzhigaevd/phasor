function saveVec2Vtk(filename,vectorStart,vectorEnd)
%  adapted by DDzhigaev for PHASOR.
%  savesthe vector to vtk

    disp('Writing Vector to VTK file....')
   
    fid = fopen(filename, 'wt');
    fprintf(fid, '# vtk DataFile Version 4.2\n');
    fprintf(fid, 'vtk output\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid, 'POINTS 1 double\n');
    fprintf(fid, '%.8f %.8f %.8f\n', vectorStart(1),vectorStart(2),vectorStart(3));
    fprintf(fid, 'CELL_TYPES 0\n');
    
    fprintf(fid, '\n');
    fprintf(fid, 'POINT_DATA 1\n');
    fprintf(fid, 'VECTORS qVector double\n');
    fprintf(fid, '%.8f %.8f %.8f\n', vectorEnd(1),vectorEnd(2),vectorEnd(3));
    fprintf(fid, '\n');
        
    fclose(fid);
    
return

