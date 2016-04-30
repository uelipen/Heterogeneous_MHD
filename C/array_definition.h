#ifndef ARRAY_DEFINITION_H
#define ARRAY_DEFINITION_H 

//      C array
/*
#define matrix4D(nub,nx,ny,nz,ii,i,j,k) ((((ii*nx)+i)*ny+j)*nz+k)
#define matrix3D(nx,ny,nz,i,j,k) (((i*ny+j)*nz)+k)
#define matrix2D(nub,nx,ii,i) (ii*nx+i)
*/

//      fortran array
#define matrix4D(nub,nx,ny,nz,ii,i,j,k) (((k*ny+j)*nx+i)*nub+ii)
#define matrix3D(nx,ny,nz,i,j,k) ((k*ny+j)*nx+i)
#define matrix2D(nub,nx,ii,i) (i*nub+ii)

#endif
