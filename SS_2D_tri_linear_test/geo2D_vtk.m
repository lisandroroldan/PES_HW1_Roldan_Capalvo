%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Getting geometry from abaqus to export it to ensight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% LOADING FILES OF NODES AND ELEM %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% WRITING HEADING FOR VTK FILE


nnode = (length(X));
nelem = length(T);

	% printing heading to file
f=fopen('MyParaviewFile.vtk','w');
fprintf(f,'# vtk DataFile Version 1.0\n');
fprintf(f,'ECM-CELL DIFFUSION-MECHANICS\n');
fprintf(f,'ASCII\n');
fprintf(f,'\n');
fprintf(f,'DATASET UNSTRUCTURED_GRID\n');
fprintf(f,'%s %8i %s\n','POINTS', nnode ,'float');

%%%%%%%%%%%%%%%%%%%%%% WRITING COORDINATES OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

	% printing coordinates
for i1=1:nnode
           fprintf(f,'%14.8E %14.8E %14.8E\n',X(i1,1:2),0.00000000E+00);
end
fprintf(f,'\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING CONNECTIVITY OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*4);
for i1=1:nelem
            fprintf(f,'%4i  %10i  %10i  %10i\n',3,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1);
end
fprintf(f,'\n');
fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
for i1=1:nelem
            fprintf(f,' %4i ', 5);
end
fprintf(f,'\n');          

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(f,'%s %8i\n','POINT_DATA', nnode);
fprintf(f,'SCALARS Ux float\n');
fprintf(f,'LOOKUP_TABLE default\n');
for i1=1:nnode
           fprintf(f,'%14.8E\n', Temp(i1) );
end

% fprintf(f,'%s %8i\n','CELL_DATA', nelem);
% fprintf(f,'SCALARS C10 float\n');
% fprintf(f,'LOOKUP_TABLE default\n');
% % for i1=1:nelem
%             fprintf(f,'%14.8E\n', Temp(i1) );
% end
fclose(f)
