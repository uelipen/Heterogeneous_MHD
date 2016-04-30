#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../data_type.h"
#include "../parameter.h"

#include "array_definition.h"
#include "value_spu.h"

#include "funclist_spu.h"

data_type_t max_3num(data_type_t *m3_1, data_type_t *m3_2, data_type_t *m3_3);
data_type_t max_2num(data_type_t *m2_1, data_type_t *m2_2);

data_type_t calcflSPE(data_type_t *calcflSPE_u, data_type_t *calcflSPE_b, int *calcflSPE_nx, int *calcflSPE_ny, int *calcflSPE_nz, int *calcflSPE_Y_location, int *calcflSPE_Z_location, int calcflSPE_tag_id)
{

volatile data_type_t calcflSPE_s_u[5*(*calcflSPE_nx)] __attribute__ ((aligned (16)));
volatile data_type_t calcflSPE_s_b[3*(*calcflSPE_nx)] __attribute__ ((aligned (16)));
volatile data_type_t calcflSPE_s_b_jp[3*(*calcflSPE_nx)] __attribute__ ((aligned (16)));
volatile data_type_t calcflSPE_s_b_kp[3*(*calcflSPE_nx)] __attribute__ ((aligned (16)));

int i,j,k; 
int i_global,j_global,k_global;

data_type_t gamma; 
gamma=5.0E0/3.0E0;
data_type_t c; 
c=1.0E0;

int ip,jp,kp; 
data_type_t bx,by,bz,v,ps,p; 
data_type_t temp1,temp2,temp3; 

int Y_inter=box_ny/num_SPE_Y;
int Z_inter=box_nz/num_SPE_Z;

int CELL_PER_BLOCK; 
CELL_PER_BLOCK = (*calcflSPE_nx);
int u_memSize_PER_BLOCK; 
u_memSize_PER_BLOCK = 5* CELL_PER_BLOCK;
int b_memSize_PER_BLOCK; 
b_memSize_PER_BLOCK = 3* CELL_PER_BLOCK;

data_type_t calcflSPE_ct;

for (k=0;k<(*calcflSPE_nz);k++)
{
        for (j=0;j<(*calcflSPE_ny);j++)
        {
//	assign global index
		k_global=k+(*calcflSPE_Z_location)*Z_inter;
		j_global=j+(*calcflSPE_Y_location)*Y_inter;
//      kp=mod(k,nz)+1
	        kp=(k_global+1)%(box_nz);
//      jp=mod(j,ny)+1
                jp=(j_global+1)%(box_ny);
//      get u
                spu_mfcdma32((void *)(calcflSPE_s_u),
                (unsigned int)(calcflSPE_u+matrix4D(5,(box_nx),(box_ny),(box_nz),0,0,(j_global),(k_global))),
                u_memSize_PER_BLOCK * sizeof(data_type_t), calcflSPE_tag_id, MFC_GETB_CMD);
//	get b
                spu_mfcdma32((void *)(calcflSPE_s_b),
                (unsigned int)(calcflSPE_b+matrix4D(3,(box_nx),(box_ny),(box_nz),0,0,(j_global),(k_global))),
                b_memSize_PER_BLOCK * sizeof(data_type_t), calcflSPE_tag_id, MFC_GET_CMD);
//      get b_jp
                spu_mfcdma32((void *)(calcflSPE_s_b_jp),
                (unsigned int)(calcflSPE_b+matrix4D(3,(box_nx),(box_ny),(box_nz),0,0,(jp),(k_global))),
                b_memSize_PER_BLOCK * sizeof(data_type_t), calcflSPE_tag_id, MFC_GET_CMD);
//      get b_kp
                spu_mfcdma32((void *)(calcflSPE_s_b_kp),
                (unsigned int)(calcflSPE_b+matrix4D(3,(box_nx),(box_ny),(box_nz),0,0,(j_global),(kp))),
                b_memSize_PER_BLOCK * sizeof(data_type_t), calcflSPE_tag_id, MFC_GET_CMD);
//      Wait for final DMAs to complete before terminating SPU thread.
                (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);
//
		for (i=0;i<(*calcflSPE_nx);i++)
		{
//	ip=mod(i,nx)+1
			i_global=i;
			ip=(i_global+1)%(box_nx);
//      bx=(b(1,i,j,k)+b(1,ip,j,k))/2.0d0
                	bx=(calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(1-1),i)]+calcflSPE_b[matrix2D(3,(*calcflSPE_nx),(1-1),ip)])/2.0E0;
//      by=(b(2,i,j,k)+b(2,i,jp,k))/2.0d0
	                by=(calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(2-1),i)]+calcflSPE_s_b_jp[matrix2D(3,(*calcflSPE_nx),(2-1),i)])/2.0E0;
//      bz=(b(3,i,j,k)+b(3,i,j,kp))/2.0d0
        	        bz=(calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(3-1),i)]+calcflSPE_s_b_kp[matrix2D(3,(*calcflSPE_nx),(3-1),i)])/2.0E0;
//      v=maxval(abs(u(2:4,i,j,k)/u(1,i,j,k)))
	                temp1=fabs(calcflSPE_s_u[matrix2D(5,(*calcflSPE_nx),(2-1),i)]/calcflSPE_s_u[matrix2D(5,(*calcflSPE_nx),(1-1),i)]);
        	        temp2=fabs(calcflSPE_s_u[matrix2D(5,(*calcflSPE_nx),(3-1),i)]/calcflSPE_s_u[matrix2D(5,(*calcflSPE_nx),(1-1),i)]);
                	temp3=fabs(calcflSPE_s_u[matrix2D(5,(*calcflSPE_nx),(4-1),i)]/calcflSPE_s_u[matrix2D(5,(*calcflSPE_nx),(1-1),i)]);
	                v=max_3num(&temp1,&temp2,&temp3);
//      ps=(u(5,i,j,k)-sum(u(2:4,i,j,k)**2,1)/u(1,i,j,k)/2.0d0)*(gamma-1.0d0)+(2.0d0-gamma)*sum(b(:,i,j,k)**2,1)/2.0d0
       		        ps=(calcflSPE_s_u[matrix2D(5,(*calcflSPE_nx),(5-1),i)]-(calcflSPE_s_u[matrix2D(5,(*calcflSPE_nx),(2-1),i)]*calcflSPE_s_u[matrix2D(5,(*calcflSPE_nx),(2-1),i)]+calcflSPE_s_u[matrix2D(5,(*calcflSPE_nx),(3-1),i)]*calcflSPE_s_u[matrix2D(5,(*calcflSPE_nx),(3-1),i)]+calcflSPE_s_u[matrix2D(5,(*calcflSPE_nx),(4-1),i)]*calcflSPE_s_u[matrix2D(5,(*calcflSPE_nx),(4-1),i)])/calcflSPE_s_u[matrix2D(5,(*calcflSPE_nx),(1-1),i)]/2.0E0)*(gamma-1.0E0)+(2.0E0-gamma)*(calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(1-1),i)]*calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(1-1),i)]+calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(2-1),i)]*calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(2-1),i)]+calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(3-1),i)]*calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(3-1),i)])/2.0E0;
//      p=ps-sum(b(:,i,j,k)**2)/2.0d0
                        p=ps-(calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(1-1),i)]*calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(1-1),i)]+calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(2-1),i)]*calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(2-1),i)]+calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(3-1),i)]*calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(3-1),i)])/2.0E0;
//      c=max(c,v+sqrt(abs(  (sum(b(:,i,j,k)**2)*2.0d0+gamma*p)/u(1,i,j,k))))
                        temp1=v+sqrt(fabs((calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(1-1),i)]*calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(1-1),i)]+calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(2-1),i)]*calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(2-1),i)]+calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(3-1),i)]*calcflSPE_s_b[matrix2D(3,(*calcflSPE_nx),(3-1),i)])*2.0E0+gamma*p)/calcflSPE_s_u[matrix2D(5,(*calcflSPE_nx),(1-1),i)]);
                        temp2=max_2num(&c,&temp1);
                        c=temp2;
		}
	}
}
calcflSPE_ct=1/c;
return calcflSPE_ct;
}

data_type_t max_3num(data_type_t *m3_1, data_type_t *m3_2, data_type_t *m3_3)
{
if ((*m3_1)>(*m3_2))
{
        if ((*m3_1)>(*m3_3))
        {
                return (*m3_1);
        }
        else
        {
                return (*m3_3);
        }
}
else
{
        if ((*m3_2)>(*m3_3))
        {
                return (*m3_2);
        }
        else
        {
                return (*m3_3);
        }
}
}

data_type_t max_2num(data_type_t *m2_1, data_type_t *m2_2)
{
if ((*m2_1)>(*m2_2))
{
        return (*m2_1);
}
else
{
        return (*m2_2);
}
}
                        
