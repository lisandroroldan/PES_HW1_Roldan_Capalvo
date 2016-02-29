%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Programming for E&S            %
%                  Homework 1                %
%       Lisandro Roldan & Albert Capalvo     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; home;
global diffusion

%Welcome message
disp(' ')
disp('This program solves a diffusion equation for the given domain.')
disp('No source term is considered.');
disp('Diffusion coeficient, geometry and BC imported from data files.');
disp(' ')

%Reading input file
disp('Example input name: Problem1');
iname = input('Input name :','s');
Xa = load(strcat('model/', iname, '_nodes'));
Ta = load(strcat('elements/', iname, '_elements'));
Inlet = importdata(strcat('model/', iname, '_groups'),',',1);
Outlet = importdata(strcat('model/', iname, '_groups'),',',3);
diffusion = load(strcat('elements/', iname, '_prop'));

%Detection of element type
nelnodes=size(Ta,2)-1;
ncoord=size(Xa,2)-1;

if ((ncoord==2) && (nelnodes==3))
 type=1;disp('Type 1: Triangular linear')
elseif ((ncoord==2) && (nelnodes==6))
 type=2;disp('Type 2: Triangular quadratic')
elseif ((ncoord==2) && (nelnodes==4))
 type=3;disp('Type 3: Cuadrilateral linear')
elseif ((ncoord==2) && (nelnodes==8))
 type=4;disp('Type 4: Cuadrilateral quadratic')
elseif ((ncoord==3) && (nelnodes==4))
 type=5;disp('Type 5: Tetra linear')
elseif ((ncoord==3) && (nelnodes==10))
 type=6;disp('Type 6: Tetra quadratic')
elseif ((ncoord==3) && (nelnodes==8))
 type=7;disp('Type 7: Hexa linear')
elseif ((ncoord==3) && (nelnodes==20))
 type=8;disp('Type 8: Hexa quadratic')
else
 disp('Wrong input')   
end

%Cutting the first column (global numbering)
X=Xa(:,2:size(Xa,2));
T=Ta(:,2:size(Ta,2));


% Geometry and BC's
% Matrix of nodal coordinates , conectivities and Boundary conditions
C = [Inlet.data(1,:)', ones(length(Inlet.data(1,:)'),1);
     Outlet.data(1,:)', zeros(length(Outlet.data(1,:)'),1)];

tic; % Start stopwatch

%Calculation of number of integrtion points
n=numberofintegrationpoints(ncoord,nelnodes);
%Position of gauss points in normalized coordinate system
pospg=integrationpoints(ncoord,nelnodes,n);
%Weights for gauss integration
wpg=integrationweights(nelnodes,n,ncoord);

%For each element the shape function and its derivative is calculated, and
%stored in matrix form.
N=zeros(n,nelnodes);
dNdxi=zeros(n*ncoord,nelnodes);   

%Loop in nodes
for i1=1:n
    Ne=shapefunctions(nelnodes,pospg,ncoord,i1);    
    N(i1,:)=Ne(:);
    dNdxie=shapefunctionderivs(nelnodes,ncoord,pospg,i1); 
    if ncoord==2
    for i2=1:nelnodes %Loop in integration points
        dNdxi(i1*2-1,i2) =dNdxie(i2,1);
        dNdxi(i1*2,i2)   =dNdxie(i2,2);
    end
    elseif ncoord==3
      for i2=1:nelnodes %Loop in integration points
        dNdxi(i1*3-2,i2) =dNdxie(i2,1);
        dNdxi(i1*3-1,i2)   =dNdxie(i2,2);
        dNdxi(i1*3,i2)   =dNdxie(i2,3);
      end   
    end
end

%Calculation of sparse matrix non-zero elements for the global K
ndir = size(C,1);
ndifzero=0;
for i=1:size(X,1)
    [row,col]=find(T == i);
    ndifzero=ndifzero+length(unique(reshape(T(row,:),[],1)));
    clear row, clear col
end
ndifzero=ndifzero+2*ndir;

% SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
[Ktot,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi,ncoord,ndifzero);

% Setting of variables needed for lagrange multiplier method
neq  = size(f,1);
A=spalloc(ndir,neq,ndir); %Redefined as sparse
A(:,C(:,1)) = eye(ndir);
b = C(:,2);

% SOLUTION OF THE LINEAR SYSTEM (using Lagrange multipliers)
% Entire matrix
Ktot = [Ktot A';A zeros(ndir,ndir)];
ftot = [f;b];
sol = Ktot\ftot;
Temp = sol(1:neq);
multip = sol(neq+1:end);

toc; %Read stopwatch


% Generation of the visualization file
name = input('Output name : ','s');
postprocess(X,T,ncoord,type,f,Temp,name);
disp('Open .vtk file in folder "result" using Paraview')

% Clear of variables 
clearvars Ne h i1 i2 iname name diffusion_in N dNdxi dNdxie Tin Xin ...
ans n b wpg pospg ncoord neq nelnodes i