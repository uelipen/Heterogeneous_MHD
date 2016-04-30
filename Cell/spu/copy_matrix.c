#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#include <stdio.h>
#include <stdlib.h>

#include "../data_type.h"
#include "../parameter.h"

#include "array_definition.h"
#include "value_spu.h"

void copy_matrix(data_type_t *copy_u, data_type_t *copy_b, int *copy_nx, int *copy_ny, int *copy_nz, data_type_t *copy_update_u, data_type_t *copy_update_b, int *copy_Y_location, int *copy_Z_location,int copy_tag_id)
{

volatile data_type_t spu_copy_u[5*(*copy_nx)];
volatile data_type_t spu_copy_b[3*(*copy_nx)];

int j,k;
int j_global,k_global;

int Y_inter=box_ny/num_SPE_Y;
int Z_inter=box_nz/num_SPE_Z;

int CELL_PER_BLOCK;
CELL_PER_BLOCK = (*copy_nx);
int u_memSize_PER_BLOCK;
u_memSize_PER_BLOCK = 5* CELL_PER_BLOCK;
int b_memSize_PER_BLOCK;
b_memSize_PER_BLOCK = 3* CELL_PER_BLOCK;

  for (k=0;k<(*copy_nz);k++)
  {
        for (j=0;j<(*copy_ny);j++)
        {
//	assign global index
		k_global=k+(*copy_Z_location)*Z_inter;
		j_global=j+(*copy_Y_location)*Y_inter;
//	get u
                spu_mfcdma32((void *)(spu_copy_u),
		(unsigned int)(copy_update_u+matrix4D(5,(box_nx),(box_ny),(box_nz),(1-1),(1-1),(j_global),(k_global))), 
		u_memSize_PER_BLOCK * sizeof(data_type_t), copy_tag_id, MFC_GETB_CMD);
//	get b
                spu_mfcdma32((void *)(spu_copy_b),
                (unsigned int)(copy_update_b+matrix4D(3,(box_nx),(box_ny),(box_nz),(1-1),(1-1),(j_global),(k_global))),
                b_memSize_PER_BLOCK * sizeof(data_type_t), copy_tag_id, MFC_GET_CMD);
//	Wait for the DMA to complete
                (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);
//	now you can copy
//
//	copy u
                spu_mfcdma32((void *)(spu_copy_u),
                (unsigned int)(copy_u+matrix4D(5,(box_nx),(box_ny),(box_nz),(1-1),(1-1),(j_global),(k_global))),
                u_memSize_PER_BLOCK * sizeof(data_type_t), copy_tag_id, MFC_PUTB_CMD);
//	copy b
                spu_mfcdma32((void *)(spu_copy_b),
                (unsigned int)(copy_b+matrix4D(3,(box_nx),(box_ny),(box_nz),(1-1),(1-1),(j_global),(k_global))),
                b_memSize_PER_BLOCK * sizeof(data_type_t), copy_tag_id, MFC_PUT_CMD);
        }
  }

//      Wait for final DMAs to complete before terminating SPU thread.
  (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);

return;
}

