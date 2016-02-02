function [n,w,xi,N,dNdxi]=T2D3
%====================== No. integration points =============================
%
%   Defines the number of integration points:be used for
%   each element type
%

npoints = 1; %1,3,4
n=npoints;
ncoord=2;  
nodes=4;
%
%====================== INTEGRATION POINTS ==================================
%
%   Defines positions of integration points

%Triangular element
       if (npoints == 1) 
         xi(1,1) = 1./3.;
         xi(2,1) = 1./3.;
       elseif (npoints == 3) 
         xi(1,1) = 0.6;
         xi(2,1) = 0.2;
         xi(1,2) = 0.2;
         xi(2,2) = 0.6;
         xi(1,3) = 0.2;
         xi(2,3) = 0.2;
       elseif (npoints == 4) 
         xi(1,1) = 1./3.;
         xi(2,1) = 1./3.;
         xi(1,2) = 0.6;
         xi(2,2) = 0.2;
         xi(1,3) = 0.2;
         xi(2,3) = 0.6;
         xi(1,4) = 0.2;
         xi(2,4) = 0.2;
       end
       
%Weights
        w = zeros(npoints,1);
       if (npoints == 1) 
         w(1) = 0.5;
       elseif (npoints == 3) 
         w(1) = 1./6.;
         w(2) = 1./6.;
         w(3) = 1./6.;
       elseif (npoints == 4) 
         w = [-27./96.,25./96.,25/96.,25/96.];
       end
%
%================= SHAPE FUNCTIONS ==================================
%
%        Nij: Shape functions of the Int Point i [4x4] Ni [4x1]

N=zeros(n,n);
for i1=1:n
       N(i1,1) = xi(i1,1);
       N(i1,2) = xi(i1,2);
       N(i1,3) = 1.-xi(i1,1)-xi(i1,2); 
end

%
%================= SHAPE FUNCTION DERIVATIVES ======================
%
%        Nij,r: Dev of shape functions of the Int Point i [4x8]
%        [2*i-1 2*i] => dNi,r [4x2]
dNdxi = zeros(ncoord*n,nodes);
for i1=1:n
       dNdxi(1,1) = 1.;
       dNdxi(1,2) = 0.;
       dNdxi(2,1) = 0.;
       dNdxi(2,2) = 1.;
       dNdxi(3,1) = -1.;
       dNdxi(3,2) = -1.;
end
end
%
