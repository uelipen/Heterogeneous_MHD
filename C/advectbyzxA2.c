#include <stdio.h>
#include <stdlib.h>

#include "array_definition.h"
#include "value_spu.h"

#include "funclist_spu.h"

void advectbyzxA2(float *advectbyzxA2_u, float *advectbyzxA2_b, int *advectbyzxA2_nx, int *advectbyzxA2_ny, int *advectbyzxA2_nz, float *advectbyzxA2_dt, float *advectbyzxA2_update_u, float *advectbyzxA2_update_b, float *advectbyzxA2_tmp)
{
int advA2_jm;
int i,j,k;

//      jm=mod(j+ny-2,ny)+1
//	b(1,:,jm,k)=b(1,:,jm,k)+cshift(fluxbx,-1)
for (k=0;k<(*advectbyzxA2_nz);k++)
{
        for (j=0;j<(*advectbyzxA2_ny);j++)
        {
                advA2_jm=(j+(*advectbyzxA2_ny)-1)%(*advectbyzxA2_ny);
		for (i=0;i<(*advectbyzxA2_nx);i++)
		{
			advectbyzxA2_b[matrix4D(3,(*advectbyzxA2_nx),(*advectbyzxA2_ny),(*advectbyzxA2_nz),(1-1),i,advA2_jm,k)]=advectbyzxA2_b[matrix4D(3,(*advectbyzxA2_nx),(*advectbyzxA2_ny),(*advectbyzxA2_nz),(1-1),i,advA2_jm,k)]+advectbyzxA2_tmp[matrix3D((*advectbyzxA2_nx),(*advectbyzxA2_ny),(*advectbyzxA2_nz),i,j,k)];
		}
	}
}

}

