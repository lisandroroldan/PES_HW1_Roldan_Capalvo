clear, close all, home
global diffusion  h 

disp(' ')
disp('This program solves a diffusion equation for the given domain.')
disp('No source term is considered.');
disp('Diffusion coeficient, geometry and BC imported from data files.');
disp(' ')

%Element type selection
disp('Element type')
disp('2D PROBLEM')
disp('1: Triangular linear')
disp('2: Triangular quadratic')
disp('3: Cuadrilateral linear')
disp('4: Cuadrilateral quadratic')
disp('3D PROBLEM')
disp('5: Tetra linear')
disp('6: Tetra quadratic')
disp('7: Hexa linear')
disp('8: Hexa quadratic')
type = input('Element type ');

%Diffusion coeficient imported
diffusion = load('Difussion/diffusion');

% GEOMETRY and BC's
% Matrix of nodal coordinates , conectivities and Boundary conditions
switch type
case 1
    Xa = load('Nodes/nodes_2D_tri_linear');
    X=Xa(:,2:3);
    Ta = load('Elem/elem_2D_tri_linear');
    T=Ta(:,2:4);
    Inlet = importdata('Groups/BC_2D_tri_linear',',',1);
    Outlet = importdata('Groups/BC_2D_tri_linear',',',3);
case 2
    Xa = load('Nodes/nodes_2D_tri_quad');
    X=Xa(:,2:3);
    Ta = load('Elem/elem_2D_tri_quad');
    T=Ta(:,2:7);
    Inlet = importdata('Groups/BC_2D_tri_quad',',',1);
    Outlet = importdata('Groups/BC_2D_tri_quad',',',3);
case 3
    Xa = load('Nodes/nodes_2D_quad_linear');
    X=Xa(:,2:3);
    Ta = load('Elem/elm_2D_quad_linear');
    T=Ta(:,2:5);
    Inlet = importdata('Groups/BC_2D_quad_linear',',',1);
    Outlet = importdata('Groups/BC_2D_quad_linear',',',3);
case 4
    Xa = load('Nodes/nodes_2D_quad_quad');
    X=Xa(:,2:3);
    Ta = load('Elem/elm_2D_quad_quad');
    T=Ta(:,2:9);
    Inlet = importdata('Groups/BC_2D_quad_quad',',',1);
    Outlet = importdata('Groups/BC_2D_quad_quad',',',3);    
case 5
    Xa = load('Nodes/nodes_3D_tri_lin');
    X=Xa(:,2:4);
    Ta = load('Elem/elem_3D_tri_lin');
    T=Ta(:,2:5);
    Inlet = importdata('Groups/BC_3D_tri_linear',',',1);
    Outlet = importdata('Groups/BC_3D_tri_linear',',',3);    
case 6
    Xa = load('Nodes/nodes_3D_tri_quad');
    X=Xa(:,2:4);
    Ta = load('Elem/elem_3D_tri_quad');
    T=Ta(:,2:11);
    Inlet = importdata('Groups/BC_3D_tri_quad',',',1);
    Outlet = importdata('Groups/BC_3D_tri_quad',',',3);    
case 7
    Xa = load('Nodes/nodes_3D_quad_lin');
    X=Xa(:,2:4);
    Ta = load('Elem/elem_3D_quad_lin');
    T=Ta(:,2:9);
    Inlet = importdata('Groups/BC_3D_quad_linear',',',1);
    Outlet = importdata('Groups/BC_3D_quad_linear',',',3);    
case 8
    Xa = load('Nodes/nodes_3D_quad_quad');
    X=Xa(:,2:4);
    Ta = load('Elem/elem_3D_quad_quad');
    T=Ta(:,2:21);
    Inlet = importdata('Groups/BC_3D_quad_quad',',',1);
    Outlet = importdata('Groups/BC_3D_quad_quad',',',3);     
otherwise
    disp('Error, non-existing element type!!! ')
end

C = [Inlet.data(1,:)', ones(length(Inlet.data(1,:)'),1);
     Outlet.data(1,:)', zeros(length(Outlet.data(1,:)'),1)];

tic;

% NUMERICAL INTEGRATION
%ncoord= Number of coordinates (2D=2, 3D=3)
%Number of nodes per element
switch type
case 1
    nelnodes=3;
    ncoord=2;
case 2
    nelnodes=6;
    ncoord=2;
case 3
    nelnodes=4;
    ncoord=2;
case 4
    nelnodes=8;
    ncoord=2;
case 5
    nelnodes=4;
    ncoord=3;
case 6
    nelnodes=10;
    ncoord=3;
case 7
    nelnodes=8;
    ncoord=3;
case 8
    nelnodes=20;
    ncoord=3;
end

%Calculation of number of integrtion points
n=numberofintegrationpoints(ncoord,nelnodes);
%Position of gauss points in normalized coordinate system
pospg=integrationpoints(ncoord,nelnodes,n);
%Weights for gauss integration
wpg=integrationweights(nelnodes,n,ncoord);

%For each element the shape function and its derivative is calculated, and
%stored in matrix form.
for i1=1:n
    Ne=shapefunctions(nelnodes,pospg,ncoord,i1);    %% this function didn't know the iteration step, thus evaluated pospg always at the same point
    N(i1,:)=Ne(:);
    dNdxie=shapefunctionderivs(nelnodes,ncoord,pospg,i1); %% same for this one, added i1
    if ncoord==2
    for i2=1:nelnodes
        dNdxi(i2,i1*2-1) =dNdxie(i2,1);
        dNdxi(i2,i1*2)   =dNdxie(i2,2);
    end
    elseif ncoord==3
      for i2=1:nelnodes
        dNdxi(i2,i1*3-2) =dNdxie(i2,1);
        dNdxi(i2,i1*3-1)   =dNdxie(i2,2);
        dNdxi(i2,i1*3)   =dNdxie(i2,3);
    end   
    end
    
    
end


%Transpose (just to make it work in CreateMatrix function
%N=N';          %No need to transpose it now
dNdxi=dNdxi' ;   % still need to transpose it, maybe tranpose it by components on lines 131-133


% SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
[K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi,ncoord);


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
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSTPROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = input('output name : ','s');
postprocess(X,T,ncoord,type,f,Temp,name);

disp('Open file: MyParaviewFile.vtk using Paraview')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear of variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars Ne h i1 i2 iname name diffusion_in N dNdxi dNdxie Tin Xin ans n b wpg pospg ncoord neq nelnodes


