function [ f ] = postprocess(X,T,ncoord,type,f,Temp,name)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nnode = length(X);
nelem = length(T);

% printing heading to file
f=fopen(strcat('results/', name, '.vtk'),'w');
fprintf(f,'# vtk DataFile Version 1.0\n');
fprintf(f,'ECM-CELL DIFFUSION-MECHANICS\n');
fprintf(f,'ASCII\n');
fprintf(f,'\n');
fprintf(f,'DATASET UNSTRUCTURED_GRID\n');
fprintf(f,'%s %8i %s\n','POINTS', nnode ,'float');

% WRITING COORDINATES OF NODES
if ncoord==2
for i1=1:nnode
fprintf(f,'%14.8E %14.8E %14.8E\n',X(i1,1:2),0.00000000E+00);
end
fprintf(f,'\n');
else
for i1=1:nnode
fprintf(f,'%14.8E %14.8E %14.8E\n',X(i1,1:3));
end
fprintf(f,'\n');
end

% WRITING CONNECTIVITY OF NODES
switch type
case 1
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
case 2
    fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*7);
    for i1=1:nelem
    fprintf(f,'%4i %10i  %10i  %10i %10i  %10i  %10i\n',6,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1);
    end
    fprintf(f,'\n');
    fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
    for i1=1:nelem
    fprintf(f,' %4i ', 22);
    end
    fprintf(f,'\n');        
case 3
    fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*5);
    for i1=1:nelem
    fprintf(f,'%4i  %10i  %10i %10i %10i\n',4,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1);
    end
    fprintf(f,'\n');
    fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
    for i1=1:nelem
    fprintf(f,' %4i ', 9);
    end
    fprintf(f,'\n');
case 4
    fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*9);
    for i1=1:nelem
    fprintf(f,'%4i %10i  %10i %10i %10i %10i  %10i %10i %10i\n',8,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1,T(i1,7)-1,T(i1,8)-1);
    end
    fprintf(f,'\n');
    fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
    for i1=1:nelem
    fprintf(f,' %4i ', 23);
    end
    fprintf(f,'\n');
case 5
    fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*5);
    for i1=1:nelem
    fprintf(f,'%4i  %10i  %10i %10i %10i\n',4,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1);
    end
    fprintf(f,'\n');
    fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
    for i1=1:nelem
    fprintf(f,' %4i ', 10);
    end
    fprintf(f,'\n'); 
case 6
    fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*11);
    for i1=1:nelem
    fprintf(f,'%4i  %10i  %10i %10i  %10i %10i  %10i %10i  %10i %10i %10i\n',10,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1,T(i1,7)-1,T(i1,8)-1,T(i1,9)-1,T(i1,10)-1);
    end
    fprintf(f,'\n');
    fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
    for i1=1:nelem
    fprintf(f,' %4i ', 24);
    end
    fprintf(f,'\n');
case 7
    fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*9);
    for i1=1:nelem
    fprintf(f,'%4i %10i  %10i %10i %10i %10i  %10i %10i %10i\n',8,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1,T(i1,7)-1,T(i1,8)-1);
    end
    fprintf(f,'\n');
    fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
    for i1=1:nelem
    fprintf(f,' %4i ', 12);
    end
    fprintf(f,'\n');
case 8
    fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*21);
    for i1=1:nelem
    fprintf(f,'%4i  %10i %10i  %10i %10i  %10i %10i %10i %10i  %10i %10i %10i %10i  %10i %10i %10i %10i %10i  %10i %10i %10i\n',20,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1,T(i1,7)-1,T(i1,8)-1,T(i1,9)-1,T(i1,10)-1,T(i1,11)-1,T(i1,12)-1,T(i1,13)-1,T(i1,14)-1,T(i1,15)-1,T(i1,16)-1,T(i1,17)-1,T(i1,18)-1,T(i1,19)-1,T(i1,20)-1);
    end
    fprintf(f,'\n');
    fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
    for i1=1:nelem
    fprintf(f,' %4i ', 25);
    end
    fprintf(f,'\n');
end

% WRITING VARIABLES
fprintf(f,'%s %8i\n','POINT_DATA', nnode);
fprintf(f,'SCALARS Ux float\n');
fprintf(f,'LOOKUP_TABLE default\n');
for i1=1:nnode
           fprintf(f,'%14.8E\n', Temp(i1) );
end

fclose(f);

end

