function N = shapefunctions(nelnodes,xi)
 

     N = zeros(nelnodes,1);

%    Triangular element

     if ( nelnodes == 3 ) 
       N(1) = xi(1);
       N(2) = xi(2);
       N(3) = 1.-xi(1)-xi(2);               
     elseif ( nelnodes == 6 ) 
       xi3 = 1.-xi(1)-xi(2);
       N(1) = (2.*xi(1)-1.)*xi(1);
       N(2) = (2.*xi(2)-1.)*xi(2);
       N(3) = (2.*xi3-1.)*xi3;
       N(4) = 4.*xi(1)*xi(2);
       N(5) = 4.*xi(2)*xi3;
       N(6) = 4.*xi3*xi(1);

%    Rectangular element
                  
     elseif ( nelnodes == 4 ) 
       N(1) = 0.25*(1.-xi(1))*(1.-xi(2));
       N(2) = 0.25*(1.+xi(1))*(1.-xi(2));
       N(3) = 0.25*(1.+xi(1))*(1.+xi(2));
       N(4) = 0.25*(1.-xi(1))*(1.+xi(2));
     elseif (nelnodes == 8) 
       N(1) = -0.25*(1.-xi(1))*(1.-xi(2))*(1.+xi(1)+xi(2));
       N(2) = 0.25*(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-1.);
       N(3) = 0.25*(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-1.);
       N(4) = 0.25*(1.-xi(1))*(1.+xi(2))*(xi(2)-xi(1)-1.);
       N(5) = 0.5*(1.-xi(1)*xi(1))*(1.-xi(2));
       N(6) = 0.5*(1.+xi(1))*(1.-xi(2)*xi(2));
       N(7) = 0.5*(1.-xi(1)*xi(1))*(1.+xi(2));
       N(8) = 0.5*(1.-xi(1))*(1.-xi(2)*xi(2));
     end


end
