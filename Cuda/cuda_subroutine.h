#ifndef CUDA_SUBROUTINE_H
#define CUDA_SUBROUTINE_H 

void cuda_cfl(float *cfl_u, float *cfl_b, int *cfl_nx, int *cfl_ny, int *cfl_nz, float *cfl_dt, float *cfl_c, int *h_cfl_nx, int *h_cfl_ny, int *h_cfl_nz, float *h_cfl_dt);

void cuda_fluidx(float *flu_u, float *flu_b, int *flu_nx, int *flu_ny, int *flu_nz, float *flu_dt, int *h_fluidx_nx, int *h_fluidx_ny, int *h_fluidx_nz);

void cuda_advectbyzx(float *adv_u, float *adv_b, int *adv_nx, int *adv_ny, int *adv_nz, float *adv_dt, float *adv_temp, int *h_adv_nx, int *h_adv_ny, int *h_adv_nz);

void cuda_transpose12(float *t12_ut, float *t12_bt, float *t12_u, float *t12_b, int *t12_nx, int *t12_ny, int *t12_nz, int *h_t12_nx, int *h_t12_ny, int *h_t12_nz);

void cuda_transpose13(float *t13_ut, float *t13_bt, float *t13_u, float *t13_b, int *t13_nx, int *t13_ny, int *t13_nz, int *h_t13_nx, int *h_t13_ny, int *h_t13_nz);

#endif

