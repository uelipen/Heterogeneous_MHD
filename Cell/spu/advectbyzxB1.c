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

void advectbyzxB1(data_type_t *advectbyzxB1_u, data_type_t *advectbyzxB1_b, int *advectbyzxB1_nx, int *advectbyzxB1_ny, int *advectbyzxB1_nz, data_type_t *advectbyzxB1_dt, data_type_t *advectbyzxB1_adv_tmp, int *advectbyzxB1_Y_location, int *advectbyzxB1_Z_location, int advectbyzxB1_tag_id)
{

volatile data_type_t advB1_s_u_km[5*(*advectbyzxB1_nx)] __attribute__ ((aligned (16)));
volatile data_type_t advB1_s_u[5*(*advectbyzxB1_nx)] __attribute__ ((aligned (16)));
volatile data_type_t advB1_s_b[3*(*advectbyzxB1_nx)] __attribute__ ((aligned (16)));

data_type_t fluxbx[(*advectbyzxB1_nx)] __attribute__ ((aligned (16)));
data_type_t b1x[(*advectbyzxB1_nx)] __attribute__ ((aligned (16)));
data_type_t vx[(*advectbyzxB1_nx)] __attribute__ ((aligned (16)));

data_type_t advB1_tmp1[(*advectbyzxB1_nx)] __attribute__ ((aligned (16)));
data_type_t advB1_tmp2[(*advectbyzxB1_nx)] __attribute__ ((aligned (16)));

int advB1_km; 
int i,j,k; 
int j_global,k_global;

int Y_inter=box_ny/num_SPE_Y;
int Z_inter=box_nz/num_SPE_Z;

int CELL_PER_BLOCK; 
CELL_PER_BLOCK = (*advectbyzxB1_nx);
int u_memSize_PER_BLOCK; 
u_memSize_PER_BLOCK = 5* CELL_PER_BLOCK;
int b_memSize_PER_BLOCK; 
b_memSize_PER_BLOCK = 3* CELL_PER_BLOCK;

for (k=0;k<(*advectbyzxB1_nz);k++)
{
	for (j=0;j<(*advectbyzxB1_ny);j++)
        {
//	assign global index
		k_global=k+(*advectbyzxB1_Z_location)*Z_inter;
		j_global=j+(*advectbyzxB1_Y_location)*Y_inter;
//	km=mod(k+nz-2,nz)+1
                advB1_km=(k_global+(box_nz)-1)%(box_nz);
//	get u(:,:,j,k)
		spu_mfcdma32((void *)(advB1_s_u),
                (unsigned int)(advectbyzxB1_u+matrix4D(5,(box_nx),(box_ny),(box_nz),0,0,(j_global),(k_global))),
                u_memSize_PER_BLOCK * sizeof(data_type_t), advectbyzxB1_tag_id, MFC_GETB_CMD);
//	get u(:,:,j,km)
                spu_mfcdma32((void *)(advB1_s_u_km),
                (unsigned int)(advectbyzxB1_u+matrix4D(5,(box_nx),(box_ny),(box_nz),0,0,(j_global),(advB1_km))),
                u_memSize_PER_BLOCK * sizeof(data_type_t), advectbyzxB1_tag_id, MFC_GET_CMD);
//	get b(:,:,j,k)
                spu_mfcdma32((void *)(advB1_s_b),
                (unsigned int)(advectbyzxB1_b+matrix4D(3,(box_nx),(box_ny),(box_nz),0,0,(j_global),(k_global))),
                b_memSize_PER_BLOCK * sizeof(data_type_t), advectbyzxB1_tag_id, MFC_GET_CMD);
//      Wait for final DMAs to complete before terminating SPU thread.
                (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);
//      vx=(u(2,:,j,km)+u(2,:,j,k))/(u(1,:,j,km)+u(1,:,j,k))
                for (i=0;i<(*advectbyzxB1_nx);i++)
                {
                        vx[i]=(advB1_s_u_km[matrix2D(5,(*advectbyzxB1_nx),(2-1),i)]+advB1_s_u[matrix2D(5,(*advectbyzxB1_nx),(2-1),i)])/(advB1_s_u_km[matrix2D(5,(*advectbyzxB1_nx),(1-1),i)]+advB1_s_u[matrix2D(5,(*advectbyzxB1_nx),(1-1),i)]);
		}
//      vx=(cshift(vx,-1)+cshift(vx,1)+2*vx)/4
                cshift_matrix1D(vx, (*advectbyzxB1_nx), (-1), advB1_tmp1);
                cshift_matrix1D(vx, (*advectbyzxB1_nx), (1), advB1_tmp2);
                for (i=0;i<(*advectbyzxB1_nx);i++)
                {
                        vx[i]=(advB1_tmp1[i]+advB1_tmp2[i]+2.0E0*vx[i])/4.0E0;
//      b1x=b(3,:,j,k)
                        b1x[i]=advB1_s_b[matrix2D(3,(*advectbyzxB1_nx),(3-1),i)];
		}
//      call tvdb(fluxbx,b1x,vx,nx,dt)
                tvdb(fluxbx,b1x,vx,advectbyzxB1_nx,advectbyzxB1_dt);
//      cshift(fluxbx,-1)
                cshift_matrix1D(fluxbx,(*advectbyzxB1_nx), (-1), advB1_tmp1);
//      b(3,:,j,k)=b1x
//      b(1,:,j,k)=b(1,:,j,k)-cshift(fluxbx,-1)
                for (i=0;i<(*advectbyzxB1_nx);i++)
                {
                        advB1_s_b[matrix2D(3,(*advectbyzxB1_nx),(3-1),i)]=b1x[i];
                        advB1_s_b[matrix2D(3,(*advectbyzxB1_nx),(1-1),i)]=advB1_s_b[matrix2D(3,(*advectbyzxB1_nx),(1-1),i)]-advB1_tmp1[i];
                }
//      send b(2,:,j,k) & b(1,:,j,k)
                spu_mfcdma32((void *)(advB1_s_b),
                (unsigned int)(advectbyzxB1_b+matrix4D(3,(box_nx),(box_ny),(box_nz),0,0,(j_global),(k_global))),
                b_memSize_PER_BLOCK * sizeof(data_type_t), advectbyzxB1_tag_id, MFC_PUTB_CMD);
//      send advB1_tmp1 out
                spu_mfcdma32((void *)(advB1_tmp1),
                (unsigned int)(advectbyzxB1_adv_tmp+matrix3D((box_nx),(box_ny),(box_nz),0,(j_global),(k_global))),
                CELL_PER_BLOCK * sizeof(data_type_t), advectbyzxB1_tag_id, MFC_PUT_CMD);
	}
}
//      Wait for final DMAs to complete before terminating SPU thread.
(void)spu_mfcstat(MFC_TAG_UPDATE_ALL);

return;
}

