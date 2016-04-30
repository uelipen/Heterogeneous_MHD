#ifndef FUNCLIST_SPU_H
#define FUNCLIST_SPU_H

void multiple_itself_matrix4D(float *mim4_ub, int mim4_nub, int mim4_nx, int mim4_ny, int mim4_nz, float mim4_multiple);

void multiple_itself_matrix2D(float *mim2_ub, int mim2_nub, int mim2_nx, float mim2_multiple);

void cshift_matrix1D(float *cm1_ub, int cm1_nx, int cm1_logical, float *cm1_out_ub);

void extract_matrix2Dto1D(float *em2to1_ub, int em2to1_nub, int em2to1_nx, int em2to1_ii, float *em2to1_out_ub);

void abXYc_matrix2D(float *abXYcm2D_L, float *abXYcm2D_Ra, float *abXYcm2D_Rb, float abXYcm2D_a, float abXYcm2D_b, int abXYcm2D_nub, int abXYcm2D_nx, float abXYcm2D_c);

void abXYc_matrix1D(float *abXYcm1D_L, float *abXYcm1D_Ra, float *abXYcm1D_Rb, float abXYcm1D_a, float abXYcm1D_b, int abXYcm1D_nx, float abXYcm1D_c);

void cshift_matrix2D_dim2(float *cm2d2_ub, int cm2d2_nub, int cm2d2_nx, int cm2d2_logical, float *cm2d2_out_ub);

void limiter_matrix2D(float *lm2D_ubA, float *lm2D_ubB, int lm2D_nub, int lm2D_nx, float *lm2D_out_ub);

#endif
