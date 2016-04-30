#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#include <stdio.h>
#include <stdlib.h>

#include "../data_type.h"
#include "../parameter.h"

#include "array_definition.h"
#include "value_spu.h"

#include "funclist_spu.h"

void advectbyzxB2(data_type_t *advectbyzxB2_b, int *advectbyzxB2_nx, int *advectbyzxB2_ny, int *advectbyzxB2_nz, data_type_t *advectbyzxB2_adv_tmp, int *advectbyzxB2_Y_location, int *advectbyzxB2_Z_location, int advectbyzxB2_tag_id)
{

volatile data_type_t advB2_s_tmp[(*advectbyzxB2_nx)] __attribute__ ((aligned (16)));
volatile data_type_t advB2_s_b_km[3*(*advectbyzxB2_nx)] __attribute__ ((aligned (16)));

int advB2_km; 
int i,j,k; 
int j_global,k_global;

int Y_inter=box_ny/num_SPE_Y;
int Z_inter=box_nz/num_SPE_Z;

int CELL_PER_BLOCK __attribute__ ((aligned (16)));
CELL_PER_BLOCK = (*advectbyzxB2_nx);
int u_memSize_PER_BLOCK __attribute__ ((aligned (16)));
u_memSize_PER_BLOCK = 5* CELL_PER_BLOCK;
int b_memSize_PER_BLOCK __attribute__ ((aligned (16)));
b_memSize_PER_BLOCK = 3* CELL_PER_BLOCK;

data_type_t temp_sizeof __attribute__ ((aligned (16)));

for (k=0;k<(*advectbyzxB2_nz);k++)
{
        for (j=0;j<(*advectbyzxB2_ny);j++)
        {
//      assign global index
                k_global=k+(*advectbyzxB2_Z_location)*Z_inter;
                j_global=j+(*advectbyzxB2_Y_location)*Y_inter;
//      km=mod(k+nz-2,nz)+1
                advB2_km=(k_global+(box_nz)-1)%(box_nz);
//	get b
                spu_mfcdma32((void *)(advB2_s_b_km),
                (unsigned int)(advectbyzxB2_b+matrix4D(3,(box_nx),(box_ny),(box_nz),0,0,(j_global),(advB2_km))),
                b_memSize_PER_BLOCK * sizeof(temp_sizeof), advectbyzxB2_tag_id, MFC_GETB_CMD);
//      get cshift(fluxbx,-1)
                spu_mfcdma32((void *)(advB2_s_tmp),
                (unsigned int)(advectbyzxB2_adv_tmp+matrix3D((box_nx),(box_ny),(box_nz),0,(j_global),(k_global))),
                CELL_PER_BLOCK * sizeof(temp_sizeof), advectbyzxB2_tag_id, MFC_GET_CMD);
//      Wait for final DMAs to complete before terminating SPU thread.
                (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);
//      b(1,:,j,km)=b(1,:,j,km)+cshift(fluxbx,-1)
		for (i=0;i<(*advectbyzxB2_nx);i++)
		{
			advB2_s_b_km[matrix2D(3,(*advectbyzxB2_nx),(1-1),i)]=advB2_s_b_km[matrix2D(3,(*advectbyzxB2_nx),(1-1),i)]+advB2_s_tmp[i];
		}
//	send b
                spu_mfcdma32((void *)(advB2_s_b_km),
                (unsigned int)(advectbyzxB2_b+matrix4D(3,(box_nx),(box_ny),(box_nz),0,0,(j_global),(advB2_km))),
                b_memSize_PER_BLOCK * sizeof(temp_sizeof), advectbyzxB2_tag_id, MFC_PUT_CMD);
	}
}
//      Wait for final DMAs to complete before terminating SPU thread.
(void)spu_mfcstat(MFC_TAG_UPDATE_ALL);

}


