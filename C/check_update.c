#include <stdio.h>
#include <stdlib.h>

#include "array_definition.h"
#include "value_spu.h"

void check_update(float *checkU_u, float *checkU_b, int *checkU_nx, int *checkU_ny, int *checkU_nz, float *checkU_update_u, float *checkU_update_b, int checkU_iter)
{
int i,j,k,ii;

int check_num=5;

	for (k=0;k<(*checkU_nz);k++)
	{
		for (j=0;j<(*checkU_ny);j++)
		{
			for (i=0;i<(*checkU_nx);i++)
			{
				for (ii=0;ii<5;ii++)
				{
if (checkU_u[matrix4D(5,(*checkU_nx),(*checkU_ny),(*checkU_nz),ii,i,j,k)]!=checkU_update_u[matrix4D(5,(*checkU_nx),(*checkU_ny),(*checkU_nz),ii,i,j,k)])
{
//	printf("update not right!\n");
//	printf("%d,%d,%d,%e,%e\n",i,j,k,checkU_u[matrix4D(5,(*checkU_nx),(*checkU_ny),(*checkU_nz),ii,i,j,k)],checkU_update_u[matrix4D(5,(*checkU_nx),(*checkU_ny),(*checkU_nz),ii,i,j,k)]);
}
				}
				for (ii=0;ii<3;ii++)
				{
if (checkU_b[matrix4D(3,(*checkU_nx),(*checkU_ny),(*checkU_nz),ii,i,j,k)]!=checkU_update_b[matrix4D(3,(*checkU_nx),(*checkU_ny),(*checkU_nz),ii,i,j,k)])
{
//        printf("update not right!\n");
//        printf("%d,%d,%d,%e,%e\n",i,j,k,checkU_b[matrix4D(3,(*checkU_nx),(*checkU_ny),(*checkU_nz),ii,i,j,k)],checkU_update_b[matrix4D(3,(*checkU_nx),(*checkU_ny),(*checkU_nz),ii,i,j,k)]);
}
				}
			}
		}
	}

}
