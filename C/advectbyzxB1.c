#include <stdio.h>
#include <stdlib.h>

#include "array_definition.h"
#include "value_spu.h"

#include "funclist_spu.h"

float fluxbx[box_nx];
float b1x[box_nx]; 
float vx[box_nx];

float advB1_tmp1[box_nx];
float advB1_tmp2[box_nx];


void tvdb(float *tvdb_fluxbx, float *tvdb_b1x, float *tvdb_vx, int *tvdb_n, float *tvdb_dt);

void advectbyzxB1(float *advectbyzxB1_u, float *advectbyzxB1_b, int *advectbyzxB1_nx, int *advectbyzxB1_ny, int *advectbyzxB1_nz, float *advectbyzxB1_dt, float *advectbyzxB1_update_u, float *advectbyzxB1_update_b, float *advectbyzxB1_tmp)
{

int advB1_km;
int i,j,k;

//	km=mod(k+nz-2,nz)+1
for (k=0;k<(*advectbyzxB1_nz);k++)
{
	for (j=0;j<(*advectbyzxB1_ny);j++)
	{
		advB1_km=(k+(*advectbyzxB1_nz)-1)%(*advectbyzxB1_nz);
		for (i=0;i<(*advectbyzxB1_nx);i++)
		{
//	vx=(u(2,:,j,km)+u(2,:,j,k))/(u(1,:,j,km)+u(1,:,j,k))
			vx[i]=(advectbyzxB1_update_u[matrix4D(5,(*advectbyzxB1_nx),(*advectbyzxB1_ny),(*advectbyzxB1_nz),(2-1),i,j,advB1_km)]+advectbyzxB1_update_u[matrix4D(5,(*advectbyzxB1_nx),(*advectbyzxB1_ny),(*advectbyzxB1_nz),(2-1),i,j,k)])/(advectbyzxB1_update_u[matrix4D(5,(*advectbyzxB1_nx),(*advectbyzxB1_ny),(*advectbyzxB1_nz),(1-1),i,j,advB1_km)]+advectbyzxB1_update_u[matrix4D(5,(*advectbyzxB1_nx),(*advectbyzxB1_ny),(*advectbyzxB1_nz),(1-1),i,j,k)]);
		}
//      vx=(cshift(vx,-1)+cshift(vx,1)+2*vx)/4
		cshift_matrix1D(vx, (*advectbyzxB1_nx), (-1), advB1_tmp1);
		cshift_matrix1D(vx, (*advectbyzxB1_nx), (1), advB1_tmp2);
                for (i=0;i<(*advectbyzxB1_nx);i++)
                {
                	vx[i]=(advB1_tmp1[i]+advB1_tmp2[i]+2.0*vx[i])/4.0;
//	b1x=b(3,:,j,k)
			b1x[i]=advectbyzxB1_update_b[matrix4D(3,(*advectbyzxB1_nx),(*advectbyzxB1_ny),(*advectbyzxB1_nz),(3-1),i,j,k)];
		}
//	call tvdb(fluxbx,b1x,vx,nx,dt)
		tvdb(fluxbx,b1x,vx,advectbyzxB1_nx,advectbyzxB1_dt);
//	cshift(fluxbx,-1)
		cshift_matrix1D(fluxbx,(*advectbyzxB1_nx), (-1), advB1_tmp1);
//	b(3,:,j,k)=b1x
//	b(1,:,j,k)=b(1,:,j,k)-cshift(fluxbx,-1)
		for (i=0;i<(*advectbyzxB1_nx);i++)
                {
			advectbyzxB1_b[matrix4D(3,(*advectbyzxB1_nx),(*advectbyzxB1_ny),(*advectbyzxB1_nz),(3-1),i,j,k)]=b1x[i];
			advectbyzxB1_b[matrix4D(3,(*advectbyzxB1_nx),(*advectbyzxB1_ny),(*advectbyzxB1_nz),(1-1),i,j,k)]=advectbyzxB1_b[matrix4D(3,(*advectbyzxB1_nx),(*advectbyzxB1_ny),(*advectbyzxB1_nz),(1-1),i,j,k)]-advB1_tmp1[i];
			advectbyzxB1_tmp[matrix3D((*advectbyzxB1_nx),(*advectbyzxB1_ny),(*advectbyzxB1_nz),i,j,k)]=advB1_tmp1[i];
		}

	}
}


}

