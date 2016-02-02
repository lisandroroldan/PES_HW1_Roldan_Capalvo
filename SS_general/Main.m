% This program solves a convection-diffusion problem 
% in a square domain [0,1]x[0,1] using bilinear elements.
% 

clear, close all, home

global diffusion  h 

disp(' ')
disp('This program solves a diffusion equation for the given domain.')
disp('No source term is considered.');
disp('Diffusion coeficient, geometry and BC imported from data files.');
disp(' ')

%Element type selection
disp('Element type')
disp('1: Triangular linear')
disp('2: Triangular quadratic')
disp('3: Cuadrilateral linear')
disp('4: Cuadrilateral quadratic')
type = input('Element type ');

%Diffusion coeficient imported
diffusion = load('diffusion');

% GEOMETRY
% Matrix of nodal coordinates and conectivities
switch type
case 1
    Xa = load('nodes_2D_tri_linear');
    X=Xa(:,2:3);
    Ta = load('elem_2D_tri_linear');
    T=Ta(:,2:4);
case 2
    Xa = load('nodes_2D_tri_quad');
    X=Xa(:,2:3);
    Ta = load('elem_2D_tri_quad');
    T=Ta(:,2:7);
case 3
    Xa = load('nodes_2D_quad_linear');
    X=Xa(:,2:3);
    Ta = load('elm_2D_quad_linear');
    T=Ta(:,2:5);
case 4
    Xa = load('nodes_2D_quad_quad');
    X=Xa(:,2:3);
    Ta = load('elm_2D_quad_quad');
    T=Ta(:,2:9);
otherwise
    disp('Error, non-existing element type!!! ')
end
 


% NUMERICAL INTEGRATION
%Number of coordinates (2D=2, 3D=3)
ncoord=2;
%Number of nodes per element
switch type
case 1
    nelnodes=3;
case 2
    nelnodes=6;
case 3
    nelnodes=4;
case 4
    nelnodes=8;
end

%Calculation of number of integrtion points
n=numberofintegrationpoints(ncoord,nelnodes);
%Position of gauss points in normalized coordinate system
pospg=integrationpoints(ncoord,nelnodes,n);
%Weights for gauss integration
wpg=integrationweights(nelnodes,n);

%For each element the shape function and its derivative is calculated, and
%stored in matrix form.
for i1=1:n
    Ne=shapefunctions(nelnodes,pospg);
    N(i1,:)=Ne(:);
    dNdxie=shapefunctionderivs(nelnodes,ncoord,pospg);
    for i2=1:nelnodes
        dNdxi(i2,i1*2-1) =dNdxie(i2,1);
        dNdxi(i2,i1*2)   =dNdxie(i2,2);
    end
end

%Transpose (just to make it work in CreateMatrix function
N=N';
dNdxi=dNdxi';

% SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
[K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi);

%BOUNDARY CONDITIONS
%nodesDir1 = nodes in wich u=1
%nodesDir2 = nodes in wich u=0
switch type
case 1
    nodesDir1 = [1, 14, 75, 76, 77, 78, 79, 80]';
    nodesDir0 = [4, 5, 6, 7, 8, 9, 10, 11, 12]';
case 2
    nodesDir1 = [1,  14,  75,  76,  77,  78,  79,  80, 453, 462, 696, 698, 701, 705, 843]'
    nodesDir0 = [4,   5,   6,   7,   8,   9,  10,  11,  12, 263, 465, 466, 468, 470, 474, 478, 711]'
case 3
    nodesDir1 = [1, 14, 75, 76, 77, 78, 79, 80]';
    nodesDir0 = [4, 5, 6, 7, 8, 9, 10, 11, 12]';
case 4
    nodesDir1 = [1,  14,  75,  76,  77,  78,  79,  80, 480, 501, 505, 607, 641, 647, 714]'
    nodesDir0 = [4,   5,   6,   7,   8,   9,  10,  11,  12, 283, 291, 487, 488, 489, 496, 728, 731]'
end

% for case 2 and 4 vectors were not transposed


% Boundary condition matrix
C = [nodesDir1, ones(length(nodesDir1),1);
     nodesDir0, zeros(length(nodesDir0),1)];

ndir = size(C,1);
neq  = size(f,1);
A = zeros(ndir,neq);
A(:,C(:,1)) = eye(ndir);
b = C(:,2);


% SOLUTION OF THE LINEAR SYSTEM
% Entire matrix
Ktot = [K A';A zeros(ndir,ndir)];
ftot = [f;b];

sol = Ktot\ftot;
Temp = sol(1:neq);
multip = sol(neq+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSTPROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nnode = length(X);
nelem = length(T);

% printing heading to file
f=fopen('MyParaviewFile.vtk','w');
fprintf(f,'# vtk DataFile Version 1.0\n');
fprintf(f,'ECM-CELL DIFFUSION-MECHANICS\n');
fprintf(f,'ASCII\n');
fprintf(f,'\n');
fprintf(f,'DATASET UNSTRUCTURED_GRID\n');
fprintf(f,'%s %8i %s\n','POINTS', nnode ,'float');

% WRITING COORDINATES OF NODES
for i1=1:nnode
fprintf(f,'%14.8E %14.8E %14.8E\n',X(i1,1:2),0.00000000E+00);
end
fprintf(f,'\n');

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
end

% WRITING VARIABLES
fprintf(f,'%s %8i\n','POINT_DATA', nnode);
fprintf(f,'SCALARS Ux float\n');
fprintf(f,'LOOKUP_TABLE default\n');
for i1=1:nnode
           fprintf(f,'%14.8E\n', Temp(i1) );
end

fclose(f);

disp('Open file: MyParaviewFile.vtk using Paraview')



