#ifndef SPU_FUNCTION_H
#define SPU_FUNCTION_H

#include "../data_type.h"

void fluidx(data_type_t *fluidx_u, data_type_t *fluidx_b, int *fluidx_nx, int *fluidx_ny, int *fluidx_nz, data_type_t *fluidx_dt, int *fluidx_Y_location, int *fluidx_Z_location, int fluidx_tag_id);

void advectbyzx(data_type_t *advectbyzx_u, data_type_t *advectbyzx_b, int *advectbyzx_nx, int *advectbyzx_ny, int *advectbyzx_nz, data_type_t *advectbyzx_dt, data_type_t *advectbyzx_adv_tmp, int *advectbyzx_Y_location, int *advectbyzx_Z_location, int advectbyzx_tag_id, int advectbyzx_SPE_id);

void transpose12(data_type_t *transpose12_u, data_type_t *transpose12_b, int *transpose12_nx, int *transpose12_ny, int *transpose12_nz, data_type_t *transpose12_trans_u_3D, data_type_t *transpose12_trans_b_3D, int *transpose12_Y_location, int *transpose12_Z_location, int transpose12_tag_id, int transpose12_SPE_id);

void transpose13(data_type_t *transpose13_u, data_type_t *transpose13_b, int *transpose13_nx, int *transpose13_ny, int *transpose13_nz, data_type_t *transpose13_trans_u_3D, data_type_t *transpose13_trans_b_3D, int *transpose13_Y_location, int *transpose13_Z_location, int transpose13_tag_id, int transpose13_SPE_id);

data_type_t calcfl(data_type_t *calcfl_u, data_type_t *calcfl_b, int *calcfl_nx, int *calcfl_ny, int *calcfl_nz, int *calcfl_Y_location, int *calcfl_Z_location, int calcfl_SPE_id, int calcfl_tag_id, unsigned int calcfl_addr);

void copy_matrix(data_type_t *copy_u, data_type_t *copy_b, int *copy_nx, int *copy_ny, int *copy_nz, data_type_t *copy_update_u, data_type_t *copy_update_b, int *copy_Y_location, int *copy_Z_location,int copy_tag_id);

void check_value(data_type_t *checkV_u, data_type_t *checkV_b, int *checkV_nx, int *checkV_ny, int *checkV_nz, int *checkV_Y_location, int *checkV_Z_location, int checkV_iter, int checkV_SPE_id, int checkV_tag_id);

//void barrier(int SPE_id, int hypercubeD);

#endif
