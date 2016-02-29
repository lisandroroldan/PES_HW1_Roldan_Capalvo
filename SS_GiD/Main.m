%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Programming for E&S            %
%     Homework 1 - GiD Problem type solver   %
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

%Import of Nodal coodrinates, conectivity, difusivity and Boundary conditions
disp('Input name example: redefined')
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


tic; %start stopwatch

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
N=zeros(n,nelnodes);
dNdxi=zeros(n*ncoord,nelnodes);   
for i1=1:n %Loop in nodes
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
ndir = size(C.data,1);
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
A=spalloc(ndir,neq,ndir); %redefined as sparse
A(:,C.data(:,1)) = eye(ndir);
b = C.data(:,2);


% SOLUTION OF THE LINEAR SYSTEM (using Lagrange multipliers)
% Entire matrix
Ktot = [Ktot A';A zeros(ndir,ndir)];
ftot = [f;b];
sol = Ktot\ftot;
Temp = sol(1:neq);
multip = sol(neq+1:end);

toc; %read stopwatch


% Generation of the visualization file
name = input('output name : ','s');
postprocess(X,T,ncoord,type,f,Temp,name);
disp('Open file: MyParaviewFile.vtk using Paraview')


% Clear of variables 
clearvars Ne h i1 i2 iname name diffusion_in N dNdxi dNdxie Tin Xin ...
ans n b wpg pospg ncoord neq nelnodes i

