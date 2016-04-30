#include <stdio.h>
#include <stdlib.h>

#include "array_definition.h"
#include "value_spu.h"

#include "funclist_spu.h"

float fluxbx[box_nx];
float b1x[box_nx]; 
float vx[box_nx];

float advA1_tmp1[box_nx];
float advA1_tmp2[box_nx];

void tvdb(float *tvdb_fluxbx, float *tvdb_b1x, float *tvdb_vx, int *tvdb_n, float *tvdb_dt);

void advectbyzxA1(float *advectbyzxA1_u, float *advectbyzxA1_b, int *advectbyzxA1_nx, int *advectbyzxA1_ny, int *advectbyzxA1_nz, float *advectbyzxA1_dt, float *advectbyzxA1_update_u, float *advectbyzxA1_update_b, float *advectbyzxA1_tmp)
{

int advA1_jm;
int i,j,k;


//	jm=mod(j+ny-2,ny)+1
for (k=0;k<(*advectbyzxA1_nz);k++)
{
	for (j=0;j<(*advectbyzxA1_ny);j++)
	{
		advA1_jm=(j+(*advectbyzxA1_ny)-1)%(*advectbyzxA1_ny);
		for (i=0;i<(*advectbyzxA1_nx);i++)
		{
//	vx=(u(2,:,jm,k)+u(2,:,j,k))/(u(1,:,jm,k)+u(1,:,j,k))
			vx[i]=(advectbyzxA1_update_u[matrix4D(5,(*advectbyzxA1_nx),(*advectbyzxA1_ny),(*advectbyzxA1_nz),(2-1),i,advA1_jm,k)]+advectbyzxA1_update_u[matrix4D(5,(*advectbyzxA1_nx),(*advectbyzxA1_ny),(*advectbyzxA1_nz),(2-1),i,j,k)])/(advectbyzxA1_update_u[matrix4D(5,(*advectbyzxA1_nx),(*advectbyzxA1_ny),(*advectbyzxA1_nz),(1-1),i,advA1_jm,k)]+advectbyzxA1_update_u[matrix4D(5,(*advectbyzxA1_nx),(*advectbyzxA1_ny),(*advectbyzxA1_nz),(1-1),i,j,k)]);
		}
//      vx=(cshift(vx,-1)+cshift(vx,1)+2*vx)/4
		cshift_matrix1D(vx, (*advectbyzxA1_nx), (-1), advA1_tmp1);
		cshift_matrix1D(vx, (*advectbyzxA1_nx), (1), advA1_tmp2);
                for (i=0;i<(*advectbyzxA1_nx);i++)
                {
                	vx[i]=(advA1_tmp1[i]+advA1_tmp2[i]+2.0*vx[i])/4.0;
//	b1x=b(2,:,j,k)
			b1x[i]=advectbyzxA1_update_b[matrix4D(3,(*advectbyzxA1_nx),(*advectbyzxA1_ny),(*advectbyzxA1_nz),(2-1),i,j,k)];
		}
//	call tvdb(fluxbx,b1x,vx,nx,dt)
		tvdb(fluxbx,b1x,vx,advectbyzxA1_nx,advectbyzxA1_dt);

//	cshift(fluxbx,-1)
		cshift_matrix1D(fluxbx,(*advectbyzxA1_nx), (-1), advA1_tmp1);
//	b(2,:,j,k)=b1x
//	b(1,:,j,k)=b(1,:,j,k)-cshift(fluxbx,-1)
		for (i=0;i<(*advectbyzxA1_nx);i++)
                {
			advectbyzxA1_b[matrix4D(3,(*advectbyzxA1_nx),(*advectbyzxA1_ny),(*advectbyzxA1_nz),(2-1),i,j,k)]=b1x[i];
			advectbyzxA1_b[matrix4D(3,(*advectbyzxA1_nx),(*advectbyzxA1_ny),(*advectbyzxA1_nz),(1-1),i,j,k)]=advectbyzxA1_b[matrix4D(3,(*advectbyzxA1_nx),(*advectbyzxA1_ny),(*advectbyzxA1_nz),(1-1),i,j,k)]-advA1_tmp1[i];
			advectbyzxA1_tmp[matrix3D((*advectbyzxA1_nx),(*advectbyzxA1_ny),(*advectbyzxA1_nz),i,j,k)]=advA1_tmp1[i];
		}

	}
}


}

