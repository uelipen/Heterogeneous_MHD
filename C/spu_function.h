void fluidx(float *fluidx_u, float *fluidx_b, int *fluidx_nx, int *fluidx_ny, int *fluidx_nz, float *fluidx_dt, float *fluidx_update_u, float *fluidx_update_b, int fluidx_iter);

void copy_matrix(float *copy_u, float *copy_b, int *copy_nx, int *copy_ny, int *copy_nz, float *copy_update_u, float *copy_update_b);

void advectbyzx(float *advectbyzx_u, float *advectbyzx_b, int *advectbyzx_nx, int *advectbyzx_ny, int *advectbyzx_nz, float *advectbyzx_dt, float *advectbyzx_update_u, float *advectbyzx_update_b);

void transpose12(float *transpose12_u, float *transpose12_b, float *transpose12_update_u, float *transpose12_update_b, int *transpose12_nx, int *transpose12_ny, int *transpose12_nz, int iter);

void transpose13(float *transpose13_u, float *transpose13_b, float *transpose13_update_u, float *transpose13_update_b, int *transpose13_nx, int *transpose13_ny, int *transpose13_nz);

void check_value(float *checkV_u, float *checkV_b, int *checkV_nx, int *checkV_ny, int *checkV_nz, float *checkV_update_u, float *checkV_update_b, int checkV_iter);

void check_update(float *checkU_u, float *checkU_b, int *checkU_nx, int *checkU_ny, int *checkU_nz, float *checkU_update_u, float *checkU_update_b, int checkU_iter);

float calcfl(float *calcfl_u, float *calcfl_b, int *calcfl_nx, int *calcfl_ny, int *calcfl_nz);

