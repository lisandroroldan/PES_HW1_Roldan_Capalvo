% This program solves a convection-diffusion problem 
% in a square domain [0,1]x[0,1] using bilinear elements.
% 

clear, close all, home

global diffusion  h 

disp(' ')
disp('This program solves a diffusion equation on [0,1]x[0,1]')
disp(' ')
disp('No source term is considered');

%  diffusion = input('Diffusion coefficient = ');
%  disp(' ')
%  nx = input('Number of elements in each direction = ');
%  ny = nx; h = 1/nx;

diffusion = load('diffusion');

% GEOMETRY
% Matrix of nodel coordinates and conectivities
%[X,T] = CreateMesh(0,1,0,1,nx,ny);
 Xa = load('nodes_2D_tri_linear');
 X=Xa(:,2:3);
 
 Ta = load('elem_2D_tri_linear');
 T=Ta(:,2:4);
 
%  nx = length(T)/sqrt(length(T));
%  ny = nx; h = 1/nx;

% NUMERICAL INTEGRATION

%User imput
nelnodes=3;
ncoord=2;
%
n=numberofintegrationpoints(ncoord,nelnodes);
pospg=integrationpoints(ncoord,nelnodes,n);
wpg=integrationweights(nelnodes,n);
for i1=1:n
    Ne=shapefunctions(nelnodes,pospg);
        N(i1,:)=Ne(:);
    dNdxie=shapefunctionderivs(nelnodes,ncoord,pospg);
        for i2=1:nelnodes
            dNdxi(i2,i1*2-1) =dNdxie(i2,1);
            dNdxi(i2,i1*2)   =dNdxie(i2,2);
        end
end
N=N';
dNdxi=dNdxi';
% Quadrature,Shape Functions
% [n,wpg,pospg,N,dNdxi] = C2D4 ;

% SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
[K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi);

% nodes on which solution is u=1
     nodesDir1 = [1, 14, 75, 76, 77, 78, 79, 80]';
     % nodes on which solution is u=0
     nodesDir0 = [4, 5, 6, 7, 8, 9, 10, 11, 12]';
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


% POSTPROCESS

