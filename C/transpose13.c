#include <stdio.h>
#include <stdlib.h>

#include "array_definition.h"
#include "value_spu.h"

void transpose12(float *transpose12_u, float *transpose12_b, float *transpose12_update_u, float *transpose12_update_b, int *transpose12_nx, int *transpose12_ny, int *transpose12_nz);

void transpose13(float *transpose13_u, float *transpose13_b, float *transpose13_update_u, float *transpose13_update_b, int *transpose13_nx, int *transpose13_ny, int *transpose13_nz)
{
int i,j,k,ii;

float f_u2_2D[5][(*transpose13_nz)][(*transpose13_nx)];
float f_b2_2D[3][(*transpose13_nz)][(*transpose13_nx)];

if ((*transpose13_nx)==(*transpose13_nz))
{
        for (j=0;j<(*transpose13_ny);j++)
        {
		for (k=0;k<(*transpose13_nz);k++)
		{
                	for (i=0;i<(*transpose13_nx);i++)
			{
			f_u2_2D[1-1][k][i]=transpose13_update_u[matrix4D(5,(*transpose13_nx),(*transpose13_ny),(*transpose13_nz),(1-1),i,j,k)];
			f_u2_2D[2-1][k][i]=transpose13_update_u[matrix4D(5,(*transpose13_nx),(*transpose13_ny),(*transpose13_nz),(4-1),i,j,k)];
			f_u2_2D[3-1][k][i]=transpose13_update_u[matrix4D(5,(*transpose13_nx),(*transpose13_ny),(*transpose13_nz),(3-1),i,j,k)];
			f_u2_2D[4-1][k][i]=transpose13_update_u[matrix4D(5,(*transpose13_nx),(*transpose13_ny),(*transpose13_nz),(2-1),i,j,k)];
			f_u2_2D[5-1][k][i]=transpose13_update_u[matrix4D(5,(*transpose13_nx),(*transpose13_ny),(*transpose13_nz),(5-1),i,j,k)];
			f_b2_2D[1-1][k][i]=transpose13_update_b[matrix4D(3,(*transpose13_nx),(*transpose13_ny),(*transpose13_nz),(3-1),i,j,k)];
			f_b2_2D[2-1][k][i]=transpose13_update_b[matrix4D(3,(*transpose13_nx),(*transpose13_ny),(*transpose13_nz),(2-1),i,j,k)];
			f_b2_2D[3-1][k][i]=transpose13_update_b[matrix4D(3,(*transpose13_nx),(*transpose13_ny),(*transpose13_nz),(1-1),i,j,k)];
			}
		}
	       	for (k=0;k<(*transpose13_nz);k++)
                {
		        for (i=0;i<(*transpose13_nx);i++)
		        {
                	        for (ii=0;ii<5;ii++)
                        	{
					transpose13_u[matrix4D(5,(*transpose13_nz),(*transpose13_ny),(*transpose13_nx),ii,k,j,i)]=f_u2_2D[ii][k][i];
				}
	                        for (ii=0;ii<3;ii++)
				{
					transpose13_b[matrix4D(3,(*transpose13_nz),(*transpose13_ny),(*transpose13_nx),ii,k,j,i)]=f_b2_2D[ii][k][i];
				}		
			}
		}
	}
}
else if ((*transpose13_nz)==1)
{
	printf("fk, what happen? nz=1?");
}
else if ((*transpose13_nx)==1)
{
	printf("fk, what happen? nx=1?");
}
else 
{
	printf("nz<>nx not supported\n");
}
}

