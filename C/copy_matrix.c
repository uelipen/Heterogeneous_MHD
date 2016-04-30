#include <stdio.h>
#include <stdlib.h>

#include <assert.h>

#include "array_definition.h"
#include "value_spu.h"

void copy_matrix(float *copy_u, float *copy_b, int *copy_nx, int *copy_ny, int *copy_nz, float *copy_update_u, float *copy_update_b)
{
int i,j,k,ii;

//printf("%d,%d,%d\n",*copy_nx,*copy_ny,*copy_nz);
  for (k=0;k<(*copy_nz);k++)
  {
        for (j=0;j<(*copy_ny);j++)
        {
		for (i=0;i<(*copy_nx);i++)
		{
			copy_update_u[matrix4D(5,(*copy_nx),(*copy_ny),(*copy_nz),(1-1),i,j,k)]=copy_u[matrix4D(5,(*copy_nx),(*copy_ny),(*copy_nz),(1-1),i,j,k)];
			copy_update_u[matrix4D(5,(*copy_nx),(*copy_ny),(*copy_nz),(2-1),i,j,k)]=copy_u[matrix4D(5,(*copy_nx),(*copy_ny),(*copy_nz),(2-1),i,j,k)];
			copy_update_u[matrix4D(5,(*copy_nx),(*copy_ny),(*copy_nz),(3-1),i,j,k)]=copy_u[matrix4D(5,(*copy_nx),(*copy_ny),(*copy_nz),(3-1),i,j,k)];
			copy_update_u[matrix4D(5,(*copy_nx),(*copy_ny),(*copy_nz),(4-1),i,j,k)]=copy_u[matrix4D(5,(*copy_nx),(*copy_ny),(*copy_nz),(4-1),i,j,k)];
			copy_update_u[matrix4D(5,(*copy_nx),(*copy_ny),(*copy_nz),(5-1),i,j,k)]=copy_u[matrix4D(5,(*copy_nx),(*copy_ny),(*copy_nz),(5-1),i,j,k)];
			copy_update_b[matrix4D(3,(*copy_nx),(*copy_ny),(*copy_nz),(1-1),i,j,k)]=copy_b[matrix4D(3,(*copy_nx),(*copy_ny),(*copy_nz),(1-1),i,j,k)];
			copy_update_b[matrix4D(3,(*copy_nx),(*copy_ny),(*copy_nz),(2-1),i,j,k)]=copy_b[matrix4D(3,(*copy_nx),(*copy_ny),(*copy_nz),(2-1),i,j,k)];
			copy_update_b[matrix4D(3,(*copy_nx),(*copy_ny),(*copy_nz),(3-1),i,j,k)]=copy_b[matrix4D(3,(*copy_nx),(*copy_ny),(*copy_nz),(3-1),i,j,k)];

if (isnan(copy_u[matrix4D(5,(*copy_nx),(*copy_ny),(*copy_nz),(5-1),i,j,k)]))
{
//	printf("NaN	%d,%d,%d\n",i,j,k);
}
		}
	}
  }

return;
}

