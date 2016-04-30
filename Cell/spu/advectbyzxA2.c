#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#include <stdio.h>
#include <stdlib.h>

#include "../data_type.h"
#include "../parameter.h"

#include "array_definition.h"
#include "value_spu.h"

#include "funclist_spu.h"

void advectbyzxA2(data_type_t *advectbyzxA2_b, int *advectbyzxA2_nx, int *advectbyzxA2_ny, int *advectbyzxA2_nz, data_type_t *advectbyzxA2_adv_tmp, int *advectbyzxA2_Y_location, int *advectbyzxA2_Z_location, int advectbyzxA2_tag_id)
{

volatile data_type_t advA2_s_tmp[(*advectbyzxA2_nx)] __attribute__ ((aligned (16)));
volatile data_type_t advA2_s_b_jm[3*(*advectbyzxA2_nx)] __attribute__ ((aligned (16)));

int advA2_jm; 
int i,j,k; 
int j_global,k_global;

int Y_inter=box_ny/num_SPE_Y;
int Z_inter=box_nz/num_SPE_Z;

int CELL_PER_BLOCK; 
CELL_PER_BLOCK = (*advectbyzxA2_nx);
int u_memSize_PER_BLOCK; 
u_memSize_PER_BLOCK = 5* CELL_PER_BLOCK;
int b_memSize_PER_BLOCK; 
b_memSize_PER_BLOCK = 3* CELL_PER_BLOCK;

for (k=0;k<(*advectbyzxA2_nz);k++)
{
        for (j=0;j<(*advectbyzxA2_ny);j++)
        {
//	assign global index
		k_global=k+(*advectbyzxA2_Z_location)*Z_inter;
		j_global=j+(*advectbyzxA2_Y_location)*Y_inter;
//      jm=mod(j+ny-2,ny)+1
                advA2_jm=(j_global+(box_ny)-1)%(box_ny);
//	get b
                spu_mfcdma32((void *)(advA2_s_b_jm),
                (unsigned int)(advectbyzxA2_b+matrix4D(3,(box_nx),(box_ny),(box_nz),0,0,(advA2_jm),(k_global))),
                b_memSize_PER_BLOCK * sizeof(data_type_t), advectbyzxA2_tag_id, MFC_GETB_CMD);
//      get cshift(fluxbx,-1)
                spu_mfcdma32((void *)(advA2_s_tmp),
                (unsigned int)(advectbyzxA2_adv_tmp+matrix3D((box_nx),(box_ny),(box_nz),0,(j_global),(k_global))),
                CELL_PER_BLOCK * sizeof(data_type_t), advectbyzxA2_tag_id, MFC_GET_CMD);
//      Wait for final DMAs to complete before terminating SPU thread.
                (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);
//      b(1,:,jm,k)=b(1,:,jm,k)+cshift(fluxbx,-1)
		for (i=0;i<(*advectbyzxA2_nx);i++)
		{
			advA2_s_b_jm[matrix2D(3,(*advectbyzxA2_nx),(1-1),i)]=advA2_s_b_jm[matrix2D(3,(*advectbyzxA2_nx),(1-1),i)]+advA2_s_tmp[i];
		}
//	send b
                spu_mfcdma32((void *)(advA2_s_b_jm),
                (unsigned int)(advectbyzxA2_b+matrix4D(3,(box_nx),(box_ny),(box_nz),0,0,(advA2_jm),(k_global))),
                b_memSize_PER_BLOCK * sizeof(data_type_t), advectbyzxA2_tag_id, MFC_PUT_CMD);
	}
}
//      Wait for final DMAs to complete before terminating SPU thread.
(void)spu_mfcstat(MFC_TAG_UPDATE_ALL);

}


