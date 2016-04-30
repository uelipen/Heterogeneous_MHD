#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#include <stdio.h>
#include <stdlib.h>

#include "../data_type.h"
#include "../parameter.h"

#include "array_definition.h"
#include "value_spu.h"

#include "funclist_spu.h"

void tvdb(data_type_t *tvdb_fluxbx, data_type_t *tvdb_b1x, data_type_t *tvdb_vx, int *tvdb_n, data_type_t *tvdb_dt);

void advectbyzxA1(data_type_t *advectbyzxA1_u, data_type_t *advectbyzxA1_b, int *advectbyzxA1_nx, int *advectbyzxA1_ny, int *advectbyzxA1_nz, data_type_t *advectbyzxA1_dt, data_type_t *advectbyzxA1_adv_tmp, int *advectbyzxA1_Y_location, int *advectbyzxA1_Z_location, int advectbyzxA1_tag_id)
{

volatile data_type_t advA1_s_u_jm[5*(*advectbyzxA1_nx)] __attribute__ ((aligned (16)));
volatile data_type_t advA1_s_u[5*(*advectbyzxA1_nx)] __attribute__ ((aligned (16)));
volatile data_type_t advA1_s_b[3*(*advectbyzxA1_nx)] __attribute__ ((aligned (16)));

data_type_t fluxbx[(*advectbyzxA1_nx)] __attribute__ ((aligned (16)));
data_type_t b1x[(*advectbyzxA1_nx)] __attribute__ ((aligned (16)));
data_type_t vx[(*advectbyzxA1_nx)] __attribute__ ((aligned (16)));

data_type_t advA1_tmp1[(*advectbyzxA1_nx)] __attribute__ ((aligned (16)));
data_type_t advA1_tmp2[(*advectbyzxA1_nx)] __attribute__ ((aligned (16)));

int advA1_jm; 
int i,j,k; 
int j_global,k_global;

int Y_inter=box_ny/num_SPE_Y;
int Z_inter=box_nz/num_SPE_Z;

int CELL_PER_BLOCK; 
CELL_PER_BLOCK = (*advectbyzxA1_nx);
int u_memSize_PER_BLOCK; 
u_memSize_PER_BLOCK = 5* CELL_PER_BLOCK;
int b_memSize_PER_BLOCK; 
b_memSize_PER_BLOCK = 3* CELL_PER_BLOCK;

for (k=0;k<(*advectbyzxA1_nz);k++)
{
	for (j=0;j<(*advectbyzxA1_ny);j++)
        {
//      assign global index
                k_global=k+(*advectbyzxA1_Z_location)*Z_inter;
                j_global=j+(*advectbyzxA1_Y_location)*Y_inter;
//      jm=mod(j+ny-2,ny)+1
		advA1_jm=(j_global+(box_ny)-1)%(box_ny);
//	get u(:,:,j,k)
		spu_mfcdma32((void *)(advA1_s_u),
                (unsigned int)(advectbyzxA1_u+matrix4D(5,(box_nx),(box_ny),(box_nz),0,0,(j_global),(k_global))),
                u_memSize_PER_BLOCK * sizeof(data_type_t), advectbyzxA1_tag_id, MFC_GETB_CMD);
//	get u(:,:,jm,k)
                spu_mfcdma32((void *)(advA1_s_u_jm),
                (unsigned int)(advectbyzxA1_u+matrix4D(5,(box_nx),(box_ny),(box_nz),0,0,(advA1_jm),(k_global))),
                u_memSize_PER_BLOCK * sizeof(data_type_t), advectbyzxA1_tag_id, MFC_GET_CMD);
//	get b(:,:,j,k)
                spu_mfcdma32((void *)(advA1_s_b),
                (unsigned int)(advectbyzxA1_b+matrix4D(3,(box_nx),(box_ny),(box_nz),0,0,(j_global),(k_global))),
                b_memSize_PER_BLOCK * sizeof(data_type_t), advectbyzxA1_tag_id, MFC_GET_CMD);
//      Wait for final DMAs to complete before terminating SPU thread.
                (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);
//      vx=(u(2,:,jm,k)+u(2,:,j,k))/(u(1,:,jm,k)+u(1,:,j,k))
		for (i=0;i<(*advectbyzxA1_nx);i++)
		{
			vx[i]=(advA1_s_u_jm[matrix2D(5,(*advectbyzxA1_nx),(2-1),i)]+advA1_s_u[matrix2D(5,(*advectbyzxA1_nx),(2-1),i)])/(advA1_s_u_jm[matrix2D(5,(*advectbyzxA1_nx),(1-1),i)]+advA1_s_u[matrix2D(5,(*advectbyzxA1_nx),(1-1),i)]);
		}
//      vx=(cshift(vx,-1)+cshift(vx,1)+2*vx)/4
                cshift_matrix1D(vx, (*advectbyzxA1_nx), (-1), advA1_tmp1);
                cshift_matrix1D(vx, (*advectbyzxA1_nx), (1), advA1_tmp2);
//	vx=(cshift(vx,-1)+cshift(vx,1)+2.0d0*vx)/4.0d0
                for (i=0;i<(*advectbyzxA1_nx);i++)
                {
                        vx[i]=(advA1_tmp1[i]+advA1_tmp2[i]+2.0E0*vx[i])/4.0E0;
			b1x[i]=advA1_s_b[matrix2D(3,(*advectbyzxA1_nx),(2-1),i)];
		}
//      call tvdb(fluxbx,b1x,vx,nx,dt)
                tvdb(fluxbx,b1x,vx,advectbyzxA1_nx,advectbyzxA1_dt);
//      cshift(fluxbx,-1)
                cshift_matrix1D(fluxbx,(*advectbyzxA1_nx), (-1), advA1_tmp1);
//      b(2,:,j,k)=b1x
//      b(1,:,j,k)=b(1,:,j,k)-cshift(fluxbx,-1)
		for (i=0;i<(*advectbyzxA1_nx);i++)
		{
			advA1_s_b[matrix2D(3,(*advectbyzxA1_nx),(2-1),i)]=b1x[i];
			advA1_s_b[matrix2D(3,(*advectbyzxA1_nx),(1-1),i)]=advA1_s_b[matrix2D(3,(*advectbyzxA1_nx),(1-1),i)]-advA1_tmp1[i];
		}
//	send b(2,:,j,k) & b(1,:,j,k)
                spu_mfcdma32((void *)(advA1_s_b),
                (unsigned int)(advectbyzxA1_b+matrix4D(3,(box_nx),(box_ny),(box_nz),0,0,(j_global),(k_global))),
                b_memSize_PER_BLOCK * sizeof(data_type_t), advectbyzxA1_tag_id, MFC_PUTB_CMD);
//	send advA1_tmp1 out
                spu_mfcdma32((void *)(advA1_tmp1),
                (unsigned int)(advectbyzxA1_adv_tmp+matrix3D((box_nx),(box_ny),(box_nz),0,(j_global),(k_global))),
                CELL_PER_BLOCK * sizeof(data_type_t), advectbyzxA1_tag_id, MFC_PUT_CMD);
	}
}
//      Wait for final DMAs to complete before terminating SPU thread.
                (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);
return;
}

