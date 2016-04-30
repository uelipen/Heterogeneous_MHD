#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#include <stdio.h>
#include <stdlib.h>

#include "../data_type.h"
#include "../parameter.h"

#include "array_definition.h"
#include "value_spu.h"

#include "funclist_spu.h"
#include "H_dma_list_elem.h"

struct dma_list_elem u_list[size_Y_trans_3D*size_Z_trans_3D] __attribute__ ((aligned (16)));
struct dma_list_elem b_list[size_Y_trans_3D*size_Z_trans_3D] __attribute__ ((aligned (16)));

void transpose13(data_type_t *transpose13_u, data_type_t *transpose13_b, int *transpose13_nx, int *transpose13_ny, int *transpose13_nz, data_type_t *transpose13_u_update, data_type_t *transpose13_b_update, int *transpose13_Y_location, int *transpose13_Z_location, int transpose13_tag_id, int transpose13_SPE_id)
{

volatile data_type_t trans_u[5*size_trans_3D] __attribute__ ((aligned (16)));
volatile data_type_t trans_b[3*size_trans_3D] __attribute__ ((aligned (16)));

volatile data_type_t trans_u_new[5*size_trans_3D] __attribute__ ((aligned (16)));
volatile data_type_t trans_b_new[3*size_trans_3D] __attribute__ ((aligned (16)));

int i,j,k;
int ii,jj,kk;

int i_global_1,j_global_1,k_global_1;
int i_global_2,j_global_2,k_global_2;

unsigned int listsize; 

int list_X=(*transpose13_nx)/size_X_trans_3D;
int list_Y=(*transpose13_ny)/size_Y_trans_3D;
int list_Z=(*transpose13_nz)/size_Z_trans_3D;

int Y_inter=box_ny/num_SPE_Y;
int Z_inter=box_nz/num_SPE_Z;

unsigned int t13_tag_id;
t13_tag_id = mfc_tag_reserve();

/*
if ((*transpose13_nx)==(*transpose13_nz))
{
*/
for (j=0;j<list_Y;j++)
{
	for (k=0;k<list_Z;k++)
	{
		for (i=0;i<list_X;i++)
		{
//	assign global index for old matrix
			i_global_1=i*size_X_trans_3D;
			j_global_1=j*size_Y_trans_3D+(*transpose13_Y_location)*Y_inter;
			k_global_1=k*size_Z_trans_3D+(*transpose13_Z_location)*Z_inter;
//	use DMA list to copy the small 3D matrix to SPE
//	set the DMA list
			for (kk=0;kk<size_Z_trans_3D;kk++)
			{
				for (jj=0;jj<size_Y_trans_3D;jj++)
				{
//
			u_list[jj+kk*size_Z_trans_3D].size.all32=5*size_X_trans_3D*sizeof(data_type_t);
			b_list[jj+kk*size_Z_trans_3D].size.all32=3*size_X_trans_3D*sizeof(data_type_t);
			u_list[jj+kk*size_Z_trans_3D].ea_low=transpose13_u+matrix4D(5,(box_nx),(box_ny),(box_nz),(1-1),(i_global_1),(j_global_1+jj),(k_global_1+kk));
			b_list[jj+kk*size_Z_trans_3D].ea_low=transpose13_b+matrix4D(3,(box_nx),(box_ny),(box_nz),(1-1),(i_global_1),(j_global_1+jj),(k_global_1+kk));
				}
			}
//      use DMA list to send u and b to SPE
                listsize=size_Y_trans_3D*size_Z_trans_3D*sizeof(struct dma_list_elem);
//      send u
                spu_mfcdma32((volatile *)(trans_u),
                (unsigned int) &u_list[0], listsize, t13_tag_id, MFC_GETLB_CMD);
//      send b
                spu_mfcdma32((volatile *)(trans_b),
                (unsigned int) &b_list[0], listsize, t13_tag_id, MFC_GETLB_CMD);
//      Wait for final DMAs to complete before terminating SPU thread.
                (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);
//	finished DMA transfer to SPE
//
//	start matrix transpose of 3D small block in SPE
			for (kk=0;kk<size_Z_trans_3D;kk++)
			{
				for (jj=0;jj<size_Y_trans_3D;jj++)
				{
					for (ii=0;ii<size_X_trans_3D;ii++)
					{
			trans_u_new[matrix4D(5,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(1-1),kk,jj,ii)]=trans_u[matrix4D(5,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(1-1),ii,jj,kk)];
			trans_u_new[matrix4D(5,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(2-1),kk,jj,ii)]=trans_u[matrix4D(5,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(4-1),ii,jj,kk)];
			trans_u_new[matrix4D(5,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(3-1),kk,jj,ii)]=trans_u[matrix4D(5,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(3-1),ii,jj,kk)];
			trans_u_new[matrix4D(5,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(4-1),kk,jj,ii)]=trans_u[matrix4D(5,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(2-1),ii,jj,kk)];
			trans_u_new[matrix4D(5,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(5-1),kk,jj,ii)]=trans_u[matrix4D(5,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(5-1),ii,jj,kk)];
			trans_b_new[matrix4D(3,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(1-1),kk,jj,ii)]=trans_b[matrix4D(3,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(3-1),ii,jj,kk)];
			trans_b_new[matrix4D(3,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(2-1),kk,jj,ii)]=trans_b[matrix4D(3,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(2-1),ii,jj,kk)];
			trans_b_new[matrix4D(3,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(3-1),kk,jj,ii)]=trans_b[matrix4D(3,size_X_trans_3D,size_Y_trans_3D,size_Z_trans_3D,(1-1),ii,jj,kk)];

					}
				}
			}
//	finished matrix transpose
//
//	assign the global index for new matrix
			i_global_2=k_global_1/size_Z_trans_3D*size_X_trans_3D;
			j_global_2=j_global_1;
			k_global_2=i_global_1/size_X_trans_3D*size_Z_trans_3D;
//	use DMA list to send u and b back to PPE
//	set the DMA list
                        for (kk=0;kk<size_Z_trans_3D;kk++)
                        {
                                for (jj=0;jj<size_Y_trans_3D;jj++)
                                {
                        u_list[jj+kk*size_Z_trans_3D].size.all32=5*size_X_trans_3D*sizeof(data_type_t);
                        b_list[jj+kk*size_Z_trans_3D].size.all32=3*size_X_trans_3D*sizeof(data_type_t);
                        u_list[jj+kk*size_Z_trans_3D].ea_low=transpose13_u_update+matrix4D(5,(box_ny),(box_nx),(box_nz),(1-1),(i_global_2),(j_global_2+jj),(k_global_2+kk));
                        b_list[jj+kk*size_Z_trans_3D].ea_low=transpose13_b_update+matrix4D(3,(box_ny),(box_nx),(box_nz),(1-1),(i_global_2),(j_global_2+jj),(k_global_2+kk));
                                }
                        }
//      use DMA list to send u and b back to PPE
                listsize=size_Y_trans_3D*size_Z_trans_3D*sizeof(struct dma_list_elem);
//      send u
                spu_mfcdma32((volatile *)(trans_u_new),
                (unsigned int) &u_list[0], listsize, t13_tag_id, MFC_PUTLB_CMD);
//      send b
                spu_mfcdma32((volatile *)(trans_b_new),
                (unsigned int) &b_list[0], listsize, t13_tag_id, MFC_PUTLB_CMD);
//      Wait for final DMAs to complete before terminating SPU thread.
                (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);
//	finished DMA transfer to PPE
		}
	}
}
//      Wait for final DMAs to complete before terminating SPU thread.
(void)spu_mfcstat(MFC_TAG_UPDATE_ALL);

mfc_tag_release(t13_tag_id);

/*
}
*/
/*
else
{
	printf("nx must be equal to nz!\n");
	return (-1);
}
*/
return;
}

