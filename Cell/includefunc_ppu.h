#ifndef INCLUDEFUNC_PPU_H
#define INCLUDEFUNC_PPU_H

#include "data_type.h"

#include "array_definition.h"

void define_value_matrix4D(data_type_t *dv4_ub, int dv4_nub, int *dv4_nx, int *dv4_ny, int *dv4_nz, data_type_t dv4_value)
{

int ii,i,j,k;
for (k=0;k<(*dv4_nz);k++)
{
	for (j=0;j<(*dv4_ny);j++)
	{
		for (i=0;i<(*dv4_nx);i++)
		{
			for (ii=0;ii<(dv4_nub);ii++)
			{
				dv4_ub[matrix4D(dv4_nub,*dv4_nx,*dv4_ny,*dv4_nz,ii,i,j,k)]=dv4_value;
			}
		}
	}
}

return;
}

void define_value_matrix3D(data_type_t *dv3_ub, int dv3_nub, int *dv3_nx, int *dv3_ny, int *dv3_nz, data_type_t dv3_value, int dv3_index)
{

int i,j,k;
for (k=0;k<(*dv3_nz);k++)
{
	for (j=0;j<(*dv3_ny);j++)
	{
		for (i=0;i<(*dv3_nx);i++)
		{
			dv3_ub[matrix4D(dv3_nub,*dv3_nx,*dv3_ny,*dv3_nz,dv3_index,i,j,k)]=dv3_value;
		}
	}
}
return;
}

void define_value_M2M_matrix4D(data_type_t *dvM2M4_ubL, data_type_t *dvM2M4_ubR, int dvM2M4_nubL, int dvM2M4_nubR, int dvM2M4_iiL, int dvM2M4_iiR, int *dvM2M4_nx, int *dvM2M4_ny, int *dvM2M4_nz, data_type_t dvM2M4_times)
{

int i,j,k;
for (k=0;k<(*dvM2M4_nz);k++)
{
	for (j=0;j<(*dvM2M4_ny);j++)
	{
		for (i=0;i<(*dvM2M4_nx);i++)
		{
			dvM2M4_ubL[matrix4D(dvM2M4_nubL,*dvM2M4_nx,*dvM2M4_ny,*dvM2M4_nz,dvM2M4_iiL,i,j,k)]=dvM2M4_times*dvM2M4_ubR[matrix4D(dvM2M4_nubR,*dvM2M4_nx,*dvM2M4_ny,*dvM2M4_nz,dvM2M4_iiR,i,j,k)];
		}	
	}
}

return;
}

#endif
