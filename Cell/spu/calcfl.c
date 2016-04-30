#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../data_type.h"
#include "../parameter.h"

data_type_t calcflSPE(data_type_t *calcflSPE_u, data_type_t *calcflSPE_b, int *calcflSPE_nx, int *calcflSPE_ny, int *calcflSPE_nz, int *calcflSPE_Y_location, int *calcflSPE_Z_location, int calcflSPE_tag_id);
data_type_t min_2num_compare(data_type_t *m2c_1, data_type_t *m2c_2);

data_type_t calcfl(data_type_t *calcfl_u, data_type_t *calcfl_b, int *calcfl_nx, int *calcfl_ny, int *calcfl_nz, int *calcfl_Y_location, int *calcfl_Z_location, int calcfl_SPE_id, int calcfl_tag_id, unsigned int *calcfl_addr)
{
data_type_t c_t_temp __attribute__ ((aligned (16)));
c_t_temp=calcflSPE(calcfl_u,calcfl_b,calcfl_nx,calcfl_ny,calcfl_nz,calcfl_Y_location,calcfl_Z_location,calcfl_tag_id);
int i;
int array_const;
array_const=16/sizeof(data_type_t);
data_type_t cfl_temp_array[MAX_SPE_THREADS*array_const] __attribute__ ((aligned (16)));
data_type_t c_t_min __attribute__ ((aligned (16)));;
unsigned int calcfl_tag;
unsigned int inside_tag_id;
inside_tag_id = mfc_tag_reserve();
barrier(calcfl_SPE_id);
if (calcfl_SPE_id!=0)
{
//	send cfl time from slave SPE to master SPE
                spu_mfcdma32((void *)(&c_t_temp),
                calcfl_addr[MAX_SPE_THREADS]+(unsigned int)&cfl_temp_array[calcfl_SPE_id*array_const],
                sizeof(data_type_t), inside_tag_id, MFC_PUTB_CMD);
}
barrier(calcfl_SPE_id);
data_type_t temp_compare;
if (calcfl_SPE_id==0)
{
	c_t_min=c_t_temp;
	for (i=1; i<MAX_SPE_THREADS; i++)
	{
//	calculate and find the smaller one
		temp_compare=min_2num_compare(&c_t_min,&cfl_temp_array[i*array_const]);	
		c_t_min=temp_compare;
	}
	for (i=1; i<MAX_SPE_THREADS; i++)
	{
                spu_mfcdma32((void *)(&c_t_min),
                calcfl_addr[MAX_SPE_THREADS+i]+(unsigned int)&c_t_min,
                sizeof(data_type_t), inside_tag_id, MFC_PUTB_CMD);
	}
//      Wait for final DMAs to complete before terminating SPU thread.
                (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);
}
barrier(calcfl_SPE_id);
mfc_tag_release(inside_tag_id);
return c_t_min;
}

data_type_t min_2num_compare(data_type_t *m2c_1, data_type_t *m2c_2)
{
if ((*m2c_1)>(*m2c_2))
{
        return (*m2c_2);
}
else
{
        return (*m2c_1);
}
}
