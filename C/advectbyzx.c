#include <stdio.h>
#include <stdlib.h>

#include "array_definition.h"
#include "value_spu.h"

float advectbyzx_tmp[box_nx*box_ny*box_nz];

void advectbyzxA1(float *advectbyzxA1_u, float *advectbyzxA1_b, int *advectbyzxA1_nx, int *advectbyzxA1_ny, int *advectbyzxA1_nz, float *advectbyzxA1_dt, float *advectbyzxA1_update_u, float *advectbyzxA1_update_b, float *advectbyzxA1_tmp);
void advectbyzxA2(float *advectbyzxA2_u, float *advectbyzxA2_b, int *advectbyzxA2_nx, int *advectbyzxA2_ny, int *advectbyzxA2_nz, float *advectbyzxA2_dt, float *advectbyzxA2_update_u, float *advectbyzxA2_update_b, float *advectbyzxA2_tmp);
void advectbyzxB1(float *advectbyzxB1_u, float *advectbyzxB1_b, int *advectbyzxB1_nx, int *advectbyzxB1_ny, int *advectbyzxB1_nz, float *advectbyzxB1_dt, float *advectbyzxB1_update_u, float *advectbyzxB1_update_b, float *advectbyzxB1_tmp);
void advectbyzxB2(float *advectbyzxB2_u, float *advectbyzxB2_b, int *advectbyzxB2_nx, int *advectbyzxB2_ny, int *advectbyzxB2_nz, float *advectbyzxB2_dt, float *advectbyzxB2_update_u, float *advectbyzxB2_update_b, float *advectbyzxB2_tmp);

void advectbyzx(float *advectbyzx_u, float *advectbyzx_b, int *advectbyzx_nx, int *advectbyzx_ny, int *advectbyzx_nz, float *advectbyzx_dt, float *advectbyzx_update_u, float *advectbyzx_update_b)
{

advectbyzxA1(advectbyzx_u,advectbyzx_b,advectbyzx_nx,advectbyzx_ny,advectbyzx_nz,advectbyzx_dt,advectbyzx_update_u,advectbyzx_update_b,advectbyzx_tmp);
advectbyzxA2(advectbyzx_u,advectbyzx_b,advectbyzx_nx,advectbyzx_ny,advectbyzx_nz,advectbyzx_dt,advectbyzx_update_u,advectbyzx_update_b,advectbyzx_tmp);
advectbyzxB1(advectbyzx_u,advectbyzx_b,advectbyzx_nx,advectbyzx_ny,advectbyzx_nz,advectbyzx_dt,advectbyzx_update_u,advectbyzx_update_b,advectbyzx_tmp);
advectbyzxB2(advectbyzx_u,advectbyzx_b,advectbyzx_nx,advectbyzx_ny,advectbyzx_nz,advectbyzx_dt,advectbyzx_update_u,advectbyzx_update_b,advectbyzx_tmp);

return;
}

