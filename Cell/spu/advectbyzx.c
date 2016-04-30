#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#include <stdio.h>
#include <stdlib.h>

#include "../data_type.h"

#include "array_definition.h"
#include "value_spu.h"

void advectbyzxA1(data_type_t *advectbyzxA1_u, data_type_t *advectbyzxA1_b, int *advectbyzxA1_nx, int *advectbyzxA1_ny, int *advectbyzxA1_nz, data_type_t *advectbyzxA1_dt, data_type_t *advectbyzxA1_adv_tmp, int *advectbyzxA1_Y_location, int *advectbyzxA1_Z_location, int advectbyzxA1_tag_id);
void advectbyzxA2(data_type_t *advectbyzxA2_b, int *advectbyzxA2_nx, int *advectbyzxA2_ny, int *advectbyzxA2_nz, data_type_t *advectbyzxA2_adv_tmp, int *advectbyzxA2_Y_location, int *advectbyzxA2_Z_location, int advectbyzxA2_tag_id);
void advectbyzxB1(data_type_t *advectbyzxB1_u, data_type_t *advectbyzxB1_b, int *advectbyzxB1_nx, int *advectbyzxB1_ny, int *advectbyzxB1_nz, data_type_t *advectbyzxB1_dt, data_type_t *advectbyzxB1_adv_tmp, int *advectbyzxB1_Y_location, int *advectbyzxB1_Z_location, int advectbyzx_tagB1_id);
void advectbyzxB2(data_type_t *advectbyzxB2_b, int *advectbyzxB2_nx, int *advectbyzxB2_ny, int *advectbyzxB2_nz, data_type_t *advectbyzxB2_adv_tmp, int *advectbyzxB2_Y_location, int *advectbyzxB2_Z_location, int advectbyzxB2_tag_id);

void advectbyzx(data_type_t *advectbyzx_u, data_type_t *advectbyzx_b, int *advectbyzx_nx, int *advectbyzx_ny, int *advectbyzx_nz, data_type_t *advectbyzx_dt, data_type_t *advectbyzx_adv_tmp, int *advectbyzx_Y_location, int *advectbyzx_Z_location, int advectbyzx_tag_id, int advectbyzx_SPE_id)
{

barrier(advectbyzx_SPE_id);
advectbyzxA1(advectbyzx_u,advectbyzx_b,advectbyzx_nx,advectbyzx_ny,advectbyzx_nz,advectbyzx_dt,advectbyzx_adv_tmp,advectbyzx_Y_location,advectbyzx_Z_location,advectbyzx_tag_id);
barrier(advectbyzx_SPE_id);
advectbyzxA2(advectbyzx_b,advectbyzx_nx,advectbyzx_ny,advectbyzx_nz,advectbyzx_adv_tmp,advectbyzx_Y_location,advectbyzx_Z_location,advectbyzx_tag_id);
barrier(advectbyzx_SPE_id);
advectbyzxB1(advectbyzx_u,advectbyzx_b,advectbyzx_nx,advectbyzx_ny,advectbyzx_nz,advectbyzx_dt,advectbyzx_adv_tmp,advectbyzx_Y_location,advectbyzx_Z_location,advectbyzx_tag_id);
barrier(advectbyzx_SPE_id);
advectbyzxB2(advectbyzx_b,advectbyzx_nx,advectbyzx_ny,advectbyzx_nz,advectbyzx_adv_tmp,advectbyzx_Y_location,advectbyzx_Z_location,advectbyzx_tag_id);
barrier(advectbyzx_SPE_id);

return;
}
