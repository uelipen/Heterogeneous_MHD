#include <stdio.h>
#include <stdlib.h>

#include "array_definition.h"
#include "value_spu.h"
#include "funclist_spu.h"
#include "includefunc_spu.h"

float spu_u[5*box_nx];
float b3x[3*box_nx];
float b2_jp[3*box_nx];
float b3_kp[3*box_nx];

float fluidx_tmp1[box_nx];
float fluidx_tmp2[box_nx]; 

void tvd1(float *tvd1_u, float *tvd1_b, int *tvd1_n, float *tvd1_dt);

void fluidx(float *fluidx_u, float *fluidx_b, int *fluidx_nx, int *fluidx_ny, int *fluidx_nz, float *fluidx_dt, float *fluidx_update_u, float *fluidx_update_b, int fluidx_iter)
{

int i,j,k,ii;
int fluidx_jp,fluidx_kp;
//
//printf("in fluidx, dt, %e\n",(*fluidx_dt));

  for (k=0;k<(*fluidx_nz);k++)
  {
        for (j=0;j<(*fluidx_ny);j++)
        {
//	get u
 		for (i=0;i<(*fluidx_nx);i++)
		{
			for (ii=0;ii<5;ii++)
			{
				spu_u[matrix2D(5,(*fluidx_nx),ii,i)]=fluidx_update_u[matrix4D(5,(*fluidx_nx),(*fluidx_ny),(*fluidx_nz),ii,i,j,k)];
if (isnan(spu_u[matrix2D(5,(*fluidx_nx),ii,i)]))
{
//        printf("in fluidx start, %d,%d,%d",i,j,k);
}
			}
		}
//	get b3x
                for (i=0;i<(*fluidx_nx);i++)
                {
                        for (ii=0;ii<3;ii++)
                        {
                                b3x[matrix2D(3,(*fluidx_nx),ii,i)]=fluidx_update_b[matrix4D(3,(*fluidx_nx),(*fluidx_ny),(*fluidx_nz),ii,i,j,k)];
                        }
                }
//	get jp & kp
		fluidx_jp=(j+1)%(*fluidx_ny);
		fluidx_kp=(k+1)%(*fluidx_nz);
//	get b(2,:,jp,k) from update data matrix
                for (i=0;i<(*fluidx_nx);i++)
                {
                        for (ii=0;ii<3;ii++)
                        {
                                b2_jp[matrix2D(3,(*fluidx_nx),ii,i)]=fluidx_update_b[matrix4D(3,(*fluidx_nx),(*fluidx_ny),(*fluidx_nz),ii,i,fluidx_jp,k)];
                        }
                }
//	get b(3,:,j,kp) from update data matrix
                for (i=0;i<(*fluidx_nx);i++)
                {
                        for (ii=0;ii<3;ii++)
                        {
                                b3_kp[matrix2D(3,(*fluidx_nx),ii,i)]=fluidx_update_b[matrix4D(3,(*fluidx_nx),(*fluidx_ny),(*fluidx_nz),ii,i,j,fluidx_kp)];
                        }
                }
//	b3x
		multiple_itself_matrix2D(b3x,3,(*fluidx_nx),0.5);
//	get cshifrt(b3x(1,:),1)
		extract_matrix2Dto1D(b3x,3,(*fluidx_nx),(1-1),fluidx_tmp1);
	 	cshift_matrix1D(fluidx_tmp1,(*fluidx_nx),1,fluidx_tmp2);	
//	b3x
		for (i=0;i<(*fluidx_nx);i++)
		{
			b3x[matrix2D(3,(*fluidx_nx),0,i)]=b3x[matrix2D(3,(*fluidx_nx),0,i)]+fluidx_tmp2[i];
			b3x[matrix2D(3,(*fluidx_nx),1,i)]=b3x[matrix2D(3,(*fluidx_nx),1,i)]+b2_jp[matrix2D(3,(*fluidx_nx),1,i)]/2;
			b3x[matrix2D(3,(*fluidx_nx),2,i)]=b3x[matrix2D(3,(*fluidx_nx),2,i)]+b3_kp[matrix2D(3,(*fluidx_nx),2,i)]/2;
		}
//	tvd1
		tvd1(spu_u,b3x,fluidx_nx,fluidx_dt);
//	copy back
                for (i=0;i<(*fluidx_nx);i++)
                {
			for (ii=0;ii<5;ii++)
                        {
				fluidx_u[matrix4D(5,(*fluidx_nx),(*fluidx_ny),(*fluidx_nz),ii,i,j,k)]=spu_u[matrix2D(5,(*fluidx_nx),ii,i)];                                
//				fluidx_u[matrix4D(5,(*fluidx_nx),(*fluidx_ny),(*fluidx_nz),ii,i,j,k)]=fluidx_u[matrix4D(5,(*fluidx_nx),(*fluidx_ny),(*fluidx_nz),ii,i,j,k)]+0.01;                                
                        }
		}
		for (i=0;i<(*fluidx_nx);i++)
                {
//			fluidx_u[matrix4D(5,(*fluidx_nx),(*fluidx_ny),(*fluidx_nz),(5-1),i,j,k)]=fluidx_u[matrix4D(5,(*fluidx_nx),(*fluidx_ny),(*fluidx_nz),(5-1),i,j,k)]+((float)(i+1))/100.0;
		}
        }
  }


return;
}

