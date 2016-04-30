#ifndef ARRAY_DEFINITION_H
#define ARRAY_DEFINITION_H 

//	fortran in C array
#define a4D_FinC(nub,nx,ny,nz,ii,i,j,k) ((((k*ny)+j)*nx+i)*nub+ii)
#define a3D_FinC(nx,ny,nz,i,j,k) (((k*ny)+j)*nx+i)
#define a2D_FinC(nub,nx,ii,i) ((i*nub)+ii) 

#endif

