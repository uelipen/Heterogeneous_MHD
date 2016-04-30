#include <stdio.h>
#include <stdlib.h>

#include "array_definition.h"
#include "value_spu.h"

#include "funclist_spu.h"

void advectbyzxB2(float *advectbyzxB2_u, float *advectbyzxB2_b, int *advectbyzxB2_nx, int *advectbyzxB2_ny, int *advectbyzxB2_nz, float *advectbyzxB2_dt, float *advectbyzxB2_update_u, float *advectbyzxB2_update_b, float *advectbyzxB2_tmp)
{
int advB2_km;
int i,j,k;

//	km=mod(k+nz-2,nz)+1
//	b(1,:,j,km)=b(1,:,j,km)+cshift(fluxbx,-1)
for (k=0;k<(*advectbyzxB2_nz);k++)
{
        for (j=0;j<(*advectbyzxB2_ny);j++)
        {
                advB2_km=(k+(*advectbyzxB2_nz)-1)%(*advectbyzxB2_nz);
		for (i=0;i<(*advectbyzxB2_nx);i++)
		{
			advectbyzxB2_b[matrix4D(3,(*advectbyzxB2_nx),(*advectbyzxB2_ny),(*advectbyzxB2_nz),(1-1),i,j,advB2_km)]=advectbyzxB2_b[matrix4D(3,(*advectbyzxB2_nx),(*advectbyzxB2_ny),(*advectbyzxB2_nz),(1-1),i,j,advB2_km)]+advectbyzxB2_tmp[matrix3D((*advectbyzxB2_nx),(*advectbyzxB2_ny),(*advectbyzxB2_nz),i,j,k)];
		}
	}
}

}


