#include <stdio.h>
#include <stdlib.h>

#include "array_definition.h"
#include "value_spu.h"

void transpose12(float *transpose12_u, float *transpose12_b, float *transpose12_update_u, float *transpose12_update_b, int *transpose12_nx, int *transpose12_ny, int *transpose12_nz, int iter)
{
int i,j,k,ii;

float t_u2_2D[5][(*transpose12_ny)][(*transpose12_nx)];
float t_b2_2D[3][(*transpose12_ny)][(*transpose12_nx)];

  FILE *check_File_fk;

for (k=0;k<(*transpose12_nz);k++)
{
        for (j=0;j<(*transpose12_ny);j++)
        {
                for (i=0;i<(*transpose12_nx);i++)
		{
			t_u2_2D[(1-1)][j][i]=transpose12_update_u[matrix4D(5,(*transpose12_nx),(*transpose12_ny),(*transpose12_nz),(1-1),i,j,k)];
			t_u2_2D[(2-1)][j][i]=transpose12_update_u[matrix4D(5,(*transpose12_nx),(*transpose12_ny),(*transpose12_nz),(3-1),i,j,k)];
			t_u2_2D[(3-1)][j][i]=transpose12_update_u[matrix4D(5,(*transpose12_nx),(*transpose12_ny),(*transpose12_nz),(2-1),i,j,k)];
			t_u2_2D[(4-1)][j][i]=transpose12_update_u[matrix4D(5,(*transpose12_nx),(*transpose12_ny),(*transpose12_nz),(4-1),i,j,k)];
			t_u2_2D[(5-1)][j][i]=transpose12_update_u[matrix4D(5,(*transpose12_nx),(*transpose12_ny),(*transpose12_nz),(5-1),i,j,k)];
			t_b2_2D[(1-1)][j][i]=transpose12_update_b[matrix4D(3,(*transpose12_nx),(*transpose12_ny),(*transpose12_nz),(2-1),i,j,k)];
			t_b2_2D[(2-1)][j][i]=transpose12_update_b[matrix4D(3,(*transpose12_nx),(*transpose12_ny),(*transpose12_nz),(1-1),i,j,k)];
			t_b2_2D[(3-1)][j][i]=transpose12_update_b[matrix4D(3,(*transpose12_nx),(*transpose12_ny),(*transpose12_nz),(3-1),i,j,k)];
		}
	}
        for (i=0;i<(*transpose12_nx);i++)
        {
        	for (j=0;j<(*transpose12_ny);j++)
                {
                        for (ii=0;ii<5;ii++)
                        {
				transpose12_u[matrix4D(5,(*transpose12_ny),(*transpose12_nx),(*transpose12_nz),ii,j,i,k)]=t_u2_2D[ii][j][i];
			}
                        for (ii=0;ii<3;ii++)
			{
				transpose12_b[matrix4D(3,(*transpose12_ny),(*transpose12_nx),(*transpose12_nz),ii,j,i,k)]=t_b2_2D[ii][j][i];
			}		
		}
	}
}
}

