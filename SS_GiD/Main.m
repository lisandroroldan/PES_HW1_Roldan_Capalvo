clear, close all, home

global diffusion  h 

disp(' ')
disp('This program solves a diffusion equation for the given domain.')
disp('No source term is considered.');
disp('Diffusion coeficient, geometry and BC imported from data files.');
disp(' ')

%Import of Nodal coodrinates, conectivity, difusivity and Boundary conditions
iname = input('input name : ','s');

header = importdata(iname,' ',1);
Xin= importdata(iname,' ',6);
Tin= importdata(iname,' ',header.data(2)+7);
diffusion_in = importdata(iname,' ',header.data(2)+8+header.data(3));
diffusion =diffusion_in.data;
C = importdata(iname,' ',header.data(2)+8+2*header.data(3)+1);
type=header.data(1);
switch type
case 1
    X=Xin.data(:,2:3);
    T=Tin.data(:,2:4);
case 2
    X=Xin.data(:,2:3);
    T=Tin.data(:,2:7);
case 3
    X=Xin.data(:,2:3);
    T=Tin.data(:,2:5);
case 4
    X=Xin.data(:,2:3);
    T=Tin.data(:,2:9);
case 5
    X=Xin.data(:,2:4);
    T=Tin.data(:,2:5);
case 6
    X=Xin.data(:,2:4);
    T=Tin.data(:,2:11);
case 7
    X=Xin.data(:,2:4);
    T=Tin.data(:,2:9);
case 8
    X=Xin.data(:,2:4);
    T=Tin.data(:,2:21);  
otherwise
    disp('Error, non-existing element type!!! ')
end


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

% Setting of variables needed for lagrange multiplier method
ndir = size(C.data,1);
neq  = size(f,1);
A = zeros(ndir,neq);
A(:,C.data(:,1)) = eye(ndir);
b = C.data(:,2);


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

