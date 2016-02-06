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
case 5
    Xa = load('nodes_3D_tri_lin');
    X=Xa(:,2:4);
    Ta = load('elem_3D_tri_lin');
    T=Ta(:,2:5);
case 6
    Xa = load('nodes_3D_tri_quad');
    X=Xa(:,2:4);
    Ta = load('elem_3D_tri_quad');
    T=Ta(:,2:11);
case 7
    Xa = load('nodes_3D_quad_lin');
    X=Xa(:,2:4);
    Ta = load('elem_3D_quad_lin');
    T=Ta(:,2:9);
case 8
    Xa = load('nodes_3D_quad_quad');
    X=Xa(:,2:4);
    Ta = load('elem_3D_quad_quad');
    T=Ta(:,2:21);  
otherwise
    disp('Error, non-existing element type!!! ')
end



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
    Ne=shapefunctions(nelnodes,pospg,ncoord);
    N(i1,:)=Ne(:);
    dNdxie=shapefunctionderivs(nelnodes,ncoord,pospg);
    for i2=1:nelnodes
        if ncoord == 2
        
        dNdxi(i2,i1*2-1) =dNdxie(i2,1);
        dNdxi(i2,i1*2)   =dNdxie(i2,2);
        else
        dNdxi(i2,i1*3-2) =dNdxie(i2,1);
        dNdxi(i2,i1*3-1) =dNdxie(i2,2);
        dNdxi(i2,i1*3)   =dNdxie(i2,3);
        
        end
        
    end
end

%Transpose (just to make it work in CreateMatrix function
N=N';
dNdxi=dNdxi';

% SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
[K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi,ncoord);

%BOUNDARY CONDITIONS
%nodesDir1 = nodes in wich u=1
%nodesDir2 = nodes in wich u=0
switch type
case 1
nodesDir1 = [1, 14, 75, 76, 77, 78, 79, 80]';
nodesDir0 = [4, 5, 6, 7, 8, 9, 10, 11, 12]';
case 2
nodesDir1 = [1,  14,  75,  76,  77,  78,  79,  80, 453, 462, 696, 698, 701, 705, 843]';
nodesDir0 = [4,   5,   6,   7,   8,   9,  10,  11,  12, 263, 465, 466, 468, 470, 474, 478, 711]';
case 3
nodesDir1 = [1, 14, 75, 76, 77, 78, 79, 80]';
nodesDir0 = [4, 5, 6, 7, 8, 9, 10, 11, 12]';
case 4
nodesDir1 = [1,  14,  75,  76,  77,  78,  79,  80, 480, 501, 505, 607, 641, 647, 714]';
nodesDir0 = [4,   5,   6,   7,   8,   9,  10,  11,  12, 283, 291, 487, 488, 489, 496, 728, 731]';
case 5
nodesDir1 = [5,   6,  11,  12,  54,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 207, 208, 209, 210, 211, 212]';
nodesDir0 = [7,   8,   9,  10,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75, 76,  77,  78,  79, 191, 192, 193, 194, 195, 196, 197]';
case 6
nodesDir1 = [5,    6,   11,   12,   54,   98,   99,  100,  101,  102,  103,  104,  105,  106,  107,  108, 109,  110,  207,  208,  209,  210,  211,  212, 1007, 1012, 1017, 1088, 1102, 1103, 1113, 1204, 1205, 1316, 1318, 1319, 1352, 1720, 1724, 1725, 1813, 1815, 1817, 1852, 1854, 1978, 2094, 2095, 2133, 2259, 2260, 2269, 2271, 2273, 2385, 2431, 2515, 2516, 2833, 2835, 2975, 3039, 3190, 3462, 3463, 3511, 3512, 3683, 3691, 3877, 4097, 4414, 4428, 4429, 4521]';
nodesDir0 = [7,    8,    9,   10,   64,   65,   66,   67,   68,   69,   70,   71,   72,   73,   74,   75,   76,   77,   78,   79,  191,  192,  193,  194,  195,  196,  197, 1033, 1133, 1136, 1141, 1168, 1174, 1239, 1241, 1243, 1248, 1311, 1315, 1486, 1497, 1506, 1668, 1669, 1761, 1763, 2810, 2857, 2859, 2860, 3468, 3469, 3697, 3936, 3946, 3948, 4114, 4116, 4117, 4194, 4196, 4228, 4229, 4249, 4452, 4455, 4456, 4461, 4463, 4464, 4465, 4466, 4467, 4470, 4471, 4472, 4474, 4475, 4477, 4484, 4485, 4486, 4489, 4490, 4507]';
case 7
nodesDir1 = [1,   3,  24,  25,  26,  27,  28,  29, 260, 262, 283, 284, 285, 286, 287, 288, 519, 521, 542, 543, 544, 545, 546, 547]';
nodesDir0 = [5,   6,   7,   8,   9,  10,  11,  12,  13, 264, 265, 266, 267, 268, 269, 270, 271, 272, 523, 524, 525, 526, 527, 528, 529, 530, 531]';
case 8
nodesDir1 = [1,    3,   24,   25,   26,   27,   28,   29,  260,  262,  283,  284,  285,  286,  287,  288,  519,  521,  542,  543,  544,  545,  546,  547,  903,  906,  909,  911,  975,  982,  984,  986,  992,  997,  998,  999, 1000, 1001, 1449, 1456, 1457, 1459, 1464, 1466, 1772, 1773, 2070, 2073, 2075, 2119, 2121, 2123, 2129, 2130, 2131, 2132, 2414, 2418, 2422, 2424, 2606]';    
nodesDir0 = [5,    6,    7,    8,    9,   10,   11,   12,   13,  264,  265,  266,  267,  268,  269,  270,  271,  272,  523,  524,  525,  526,  527,  528,  529,  530,  531,  839,  844,  848,  849,  852,  855,  858,  860,  862,  865, 1413, 1420, 1422, 1423, 1424, 1431, 1432, 1434, 1974, 1977, 1978, 1979, 1982, 1983, 1984, 2031, 2035, 2036, 2038, 2041, 2043, 2046, 2395, 2397, 2398, 2402, 2403, 2405, 2717, 2718, 2720, 2721]';            
end


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


F= postprocess(X,T,ncoord,type,f,Temp);

disp('Open file: MyParaviewFile.vtk using Paraview')



