*loop elems
*set var el=ElemsType
*set var quad=IsQuadratic
*if(el==2 && quad==0 )
*set var type=1
*elseif(el==2 && quad==1 )
*set var type=2
*elseif(el==3 && quad==0 )
*set var type=3
*elseif(el==3 && quad==1 )
*set var type=4
*elseif(el==4 && quad==0 )
*set var type=5
*elseif(el==4 && quad==1 )
*set var type=6
*elseif(el==5 && quad==0 )
*set var type=7
*elseif(el==5 && quad==1 )
*set var type=7
*end fi
*# ///////////////////// linear triangles//////////////////////////////////////////////////////
*# 
*if(type==1)
elemt type / number nodes / number elements / number conditions
*type
*npoin
*nelem
*Set Cond Dirichlet_BC_line *nodes
*CondNumEntities
Nodal coordinate matrix /////////////////////////////////
*loop nodes
*NodesNum *NodesCoord(1) *NodesCoord(2)
*end nodes
COnectivity matrix /////////////////////////////////
*loop elems 
*ElemsNum *ElemsConec(1) *ElemsConec(2) *ElemsConec(3)
*end elems
material properties ///////////////////////////////////////
*loop materials
*loop elems 
*if(matnum==elemsmat)
*ElemsNum  *MatProp(1) 
*end fi
*end elems
*end materials
Boundary conditions //////////////////////////////////////////
*Set Cond Dirichlet_BC_line *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(1)  
*end nodes
*end fi
*# ///////////////////// quadratic triangles////////////////////////////////////////////////////////////
*# 
*if(type==2)
elemt type / number nodes / number elements / number conditions
*type
*npoin
*nelem
*Set Cond Dirichlet_BC_line *nodes
*CondNumEntities
Nodal coordinate matrix /////////////////////////////////
*loop nodes
*NodesNum *NodesCoord(1) *NodesCoord(2)
*end nodes
COnectivity matrix /////////////////////////////////
*loop elems 
*ElemsNum *ElemsConec(1) *ElemsConec(2) *ElemsConec(3) *ElemsConec(4) *ElemsConec(5) *ElemsConec(6)
*end elems
material properties ///////////////////////////////////////
*loop materials
*loop elems 
*if(matnum==elemsmat)
*ElemsNum  *MatProp(1) 
*end fi
*end elems
*end materials
Boundary conditions //////////////////////////////////////////
*Set Cond Dirichlet_BC_line *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(1)  
*end nodes
*end fi
*# 
*# ///////////////////// linear quadrilateral///////////////////////////////////////////////////////////
*# 
*if(type==3)
elemt type / number nodes / number elements / number conditions
*type
*npoin
*nelem
*Set Cond Dirichlet_BC_line *nodes
*CondNumEntities
Nodal coordinate matrix /////////////////////////////////
*loop nodes
*NodesNum *NodesCoord(1) *NodesCoord(2)
*end nodes
COnectivity matrix /////////////////////////////////
*loop elems 
*ElemsNum *ElemsConec(1) *ElemsConec(2) *ElemsConec(3) *ElemsConec(4)
*end elems
material properties ///////////////////////////////////////
*loop materials
*loop elems 
*if(matnum==elemsmat)
*ElemsNum  *MatProp(1) 
*end fi
*end elems
*end materials
Boundary conditions //////////////////////////////////////////
*Set Cond Dirichlet_BC_line *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(1)  
*end nodes
*end fi
*# 
*# ///////////////////// quadratic quadrilateral///////////////////////////////////////////////////////////
*# 
*if(type==4)
elemt type / number nodes / number elements / number conditions
*type
*npoin
*nelem
*Set Cond Dirichlet_BC_line *nodes
*CondNumEntities
Nodal coordinate matrix /////////////////////////////////
*loop nodes
*NodesNum *NodesCoord(1) *NodesCoord(2)
*end nodes
COnectivity matrix /////////////////////////////////
*loop elems 
*ElemsNum *ElemsConec(1) *ElemsConec(2) *ElemsConec(3) *ElemsConec(4) *ElemsConec(5) *ElemsConec(6) *ElemsConec(7) *ElemsConec(8)
*end elems
material properties ///////////////////////////////////////
*loop materials
*loop elems 
*if(matnum==elemsmat)
*ElemsNum  *MatProp(1) 
*end fi
*end elems
*end materials
Boundary conditions //////////////////////////////////////////
*Set Cond Dirichlet_BC_line *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(1)  
*end nodes
*end fi
*# 
*# ///////////////////// linear tetahedra///////////////////////////////////////////////////////////
*# 
*if(type==5)
elemt type / number nodes / number elements / number conditions
*type
*npoin
*nelem
*Set Cond Dirichlet_BC_line *nodes
*add cond Dirichlet_BC_surface *nodes
*CondNumEntities
Nodal coordinate matrix /////////////////////////////////
*loop nodes
*NodesNum *NodesCoord(1) *NodesCoord(2) *NodesCoord(3)
*end nodes
COnectivity matrix /////////////////////////////////
*loop elems 
*ElemsNum *ElemsConec(1) *ElemsConec(2) *ElemsConec(3) *ElemsConec(4)
*end elems
material properties ///////////////////////////////////////
*loop materials
*loop elems 
*if(matnum==elemsmat)
*ElemsNum  *MatProp(1) 
*end fi
*end elems
*end materials
Boundary conditions //////////////////////////////////////////
*Set Cond Dirichlet_BC_line *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(1)  
*end nodes
*Set Cond Dirichlet_BC_surface *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(1)  
*end nodes
*end fi
*# 
*# ///////////////////// quadratic tetahedra///////////////////////////////////////////////////////////
*# 
*if(type==6)
elemt type / number nodes / number elements / number conditions
*type
*npoin
*nelem
*Set Cond Dirichlet_BC_line *nodes
*add cond Dirichlet_BC_surface *nodes
*CondNumEntities
Nodal coordinate matrix /////////////////////////////////
*loop nodes
*NodesNum *NodesCoord(1) *NodesCoord(2) *NodesCoord(3)
*end nodes
COnectivity matrix /////////////////////////////////
*loop elems 
*ElemsNum *ElemsConec(1) *ElemsConec(2) *ElemsConec(3) *ElemsConec(4) *ElemsConec(5) *ElemsConec(6) *ElemsConec(7) *ElemsConec(8) *ElemsConec(9) *ElemsConec(10)
*end elems
material properties ///////////////////////////////////////
*loop materials
*loop elems 
*if(matnum==elemsmat)
*ElemsNum  *MatProp(1) 
*end fi
*end elems
*end materials
Boundary conditions //////////////////////////////////////////
*Set Cond Dirichlet_BC_line *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(1)  
*end nodes
*Set Cond Dirichlet_BC_surface *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(1)  
*end nodes
*end fi
*# 
*# ///////////////////// linear hexahedron///////////////////////////////////////////////////////////
*# 
*if(type==7)
elemt type / number nodes / number elements / number conditions
*type
*npoin
*nelem
*Set Cond Dirichlet_BC_line *nodes
*add cond Dirichlet_BC_surface *nodes
*CondNumEntities
Nodal coordinate matrix /////////////////////////////////
*loop nodes
*NodesNum *NodesCoord(1) *NodesCoord(2) *NodesCoord(3)
*end nodes
COnectivity matrix /////////////////////////////////
*loop elems 
*ElemsNum *ElemsConec(1) *ElemsConec(2) *ElemsConec(3) *ElemsConec(4) *ElemsConec(5) *ElemsConec(6) *ElemsConec(7) *ElemsConec(8)
*end elems
material properties ///////////////////////////////////////
*loop materials
*loop elems 
*if(matnum==elemsmat)
*ElemsNum  *MatProp(1) 
*end fi
*end elems
*end materials
Boundary conditions //////////////////////////////////////////
*Set Cond Dirichlet_BC_line *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(1)  
*end nodes
*Set Cond Dirichlet_BC_surface *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(1)  
*end nodes
*end fi
*# 
*# ///////////////////// quadratic hexahedron///////////////////////////////////////////////////////////
*# 
*if(type==8)
elemt type / number nodes / number elements / number conditions
*type
*npoin
*nelem
*Set Cond Dirichlet_BC_line *nodes
*add cond Dirichlet_BC_surface *nodes
*CondNumEntities
Nodal coordinate matrix /////////////////////////////////
*loop nodes
*NodesNum *NodesCoord(1) *NodesCoord(2) *NodesCoord(3)
*end nodes
COnectivity matrix /////////////////////////////////
*loop elems 
*ElemsNum *ElemsConec(1) *ElemsConec(2) *ElemsConec(3) *ElemsConec(4) *ElemsConec(5) *ElemsConec(6) *ElemsConec(7) *ElemsConec(8) *ElemsConec(9) *ElemsConec(10) *ElemsConec(11) *ElemsConec(12) *ElemsConec(13) *ElemsConec(14) *ElemsConec(15) *ElemsConec(16) *ElemsConec(17) *ElemsConec(18) *ElemsConec(19) *ElemsConec(20) 
*end elems
material properties ///////////////////////////////////////
*loop materials
*loop elems 
*if(matnum==elemsmat)
*ElemsNum  *MatProp(1) 
*end fi
*end elems
*end materials
Boundary conditions //////////////////////////////////////////
*Set Cond Dirichlet_BC_line *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(1)  
*end nodes
*Set Cond Dirichlet_BC_surface *nodes
*loop nodes *OnlyInCond
*NodesNum *cond(1)  
*end nodes
*end fi
