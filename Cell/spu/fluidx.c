#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#include <stdio.h>
#include <stdlib.h>

#include "../data_type.h"
#include "../parameter.h"

#include "array_definition.h"
#include "value_spu.h"
#include "funclist_spu.h"
#include "includefunc_spu.h"

void tvd1(data_type_t *tvd1_u, data_type_t *tvd1_b, int *tvd1_n, data_type_t *tvd1_dt);

void fluidx(data_type_t *fluidx_u, data_type_t *fluidx_b, int *fluidx_nx, int *fluidx_ny, int *fluidx_nz, data_type_t *fluidx_dt, int *fluidx_Y_location, int *fluidx_Z_location, int fluidx_tag_id)
{

volatile data_type_t fluidx_s_u[5*(*fluidx_nx)] __attribute__ ((aligned (16)));
volatile data_type_t fluidx_s_b3x[3*(*fluidx_nx)] __attribute__ ((aligned (16)));
volatile data_type_t fluidx_s_b2_jp[3*(*fluidx_nx)] __attribute__ ((aligned (16)));
volatile data_type_t fluidx_s_b3_kp[3*(*fluidx_nx)] __attribute__ ((aligned (16)));

data_type_t fluidx_tmp1[(*fluidx_nx)] __attribute__ ((aligned (16)));
data_type_t fluidx_tmp2[(*fluidx_nx)] __attribute__ ((aligned (16)));

int i,j,k; 
int j_global,k_global;
int fluidx_jp,fluidx_kp; 

int Y_inter=box_ny/num_SPE_Y;
int Z_inter=box_nz/num_SPE_Z;

int CELL_PER_BLOCK; 
CELL_PER_BLOCK = (*fluidx_nx);
int u_memSize_PER_BLOCK; 
u_memSize_PER_BLOCK = 5* CELL_PER_BLOCK; 
int b_memSize_PER_BLOCK; 
b_memSize_PER_BLOCK = 3* CELL_PER_BLOCK;

  for (k=0;k<(*fluidx_nz);k++)
  {
        for (j=0;j<(*fluidx_ny);j++)
        {
//	assigne global index
		k_global=k+(*fluidx_Z_location)*Z_inter;
		j_global=j+(*fluidx_Y_location)*Y_inter;
//      get u
                spu_mfcdma32((void *)(fluidx_s_u),
                (unsigned int)(fluidx_u+matrix4D(5,(box_nx),(box_ny),(box_nz),0,0,(j_global),(k_global))),
                u_memSize_PER_BLOCK * sizeof(data_type_t), fluidx_tag_id, MFC_GETB_CMD);
//      get b3x
                spu_mfcdma32((void *)(fluidx_s_b3x),
                (unsigned int)(fluidx_b+matrix4D(3,(box_nx),(box_ny),(box_nz),0,0,(j_global),(k_global))),
                b_memSize_PER_BLOCK * sizeof(data_type_t), fluidx_tag_id, MFC_GET_CMD);
//	get jp & kp
                fluidx_jp=(j_global+1)%(box_ny);
                fluidx_kp=(k_global+1)%(box_nz);
//	get b(2,:,jp,k) from update data matrix
                spu_mfcdma32((void *)(fluidx_s_b2_jp),
                (unsigned int)(fluidx_b+matrix4D(3,(box_nx),(box_ny),(box_nz),0,0,(fluidx_jp),(k_global))),
                b_memSize_PER_BLOCK * sizeof(data_type_t), fluidx_tag_id, MFC_GET_CMD);
//	get b(3,:,j,kp) from update data matrix

                spu_mfcdma32((void *)(fluidx_s_b3_kp),
                (unsigned int)(fluidx_b+matrix4D(3,(box_nx),(box_ny),(box_nz),0,0,(j_global),(fluidx_kp))),
                b_memSize_PER_BLOCK * sizeof(data_type_t), fluidx_tag_id, MFC_GET_CMD);
//      Wait for final DMAs to complete before terminating SPU thread.
		(void)spu_mfcstat(MFC_TAG_UPDATE_ALL);
//	b3x
		multiple_itself_matrix2D(fluidx_s_b3x,3,(*fluidx_nx),0.5E0);
//	get cshifrt(fluidx_s_b3x(1,:),1)
		extract_matrix2Dto1D(fluidx_s_b3x,3,(*fluidx_nx),(1-1),fluidx_tmp1);
	 	cshift_matrix1D(fluidx_tmp1,(*fluidx_nx),1,fluidx_tmp2);	
//	b3x
		for (i=0;i<(*fluidx_nx);i++)
		{
			fluidx_s_b3x[matrix2D(3,(*fluidx_nx),0,i)]=fluidx_s_b3x[matrix2D(3,(*fluidx_nx),0,i)]+fluidx_tmp2[i];
			fluidx_s_b3x[matrix2D(3,(*fluidx_nx),1,i)]=fluidx_s_b3x[matrix2D(3,(*fluidx_nx),1,i)]+fluidx_s_b2_jp[matrix2D(3,(*fluidx_nx),1,i)]/2E0;
			fluidx_s_b3x[matrix2D(3,(*fluidx_nx),2,i)]=fluidx_s_b3x[matrix2D(3,(*fluidx_nx),2,i)]+fluidx_s_b3_kp[matrix2D(3,(*fluidx_nx),2,i)]/2E0;
		}
//	tvd1
		tvd1(fluidx_s_u,fluidx_s_b3x,fluidx_nx,fluidx_dt);
//	send u
                spu_mfcdma32((void *)(fluidx_s_u),
                (unsigned int)(fluidx_u+matrix4D(5,(box_nx),(box_ny),(box_nz),0,0,(j_global),(k_global))),
                u_memSize_PER_BLOCK * sizeof(data_type_t), fluidx_tag_id, MFC_PUT_CMD);
        }
  }
//      Wait for final DMAs to complete before terminating SPU thread.
  (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);

return;
}

