#ifndef FUNCLIST_SPU_H
#define FUNCLIST_SPU_H

#include "../data_type.h"

void multiple_itself_matrix4D(data_type_t *mim4_ub, int mim4_nub, int mim4_nx, int mim4_ny, int mim4_nz, data_type_t mim4_multiple);

void multiple_itself_matrix2D(data_type_t *mim2_ub, int mim2_nub, int mim2_nx, data_type_t mim2_multiple);

void cshift_matrix1D(data_type_t *cm1_ub, int cm1_nx, int cm1_logical, data_type_t *cm1_out_ub);

void extract_matrix2Dto1D(data_type_t *em2to1_ub, int em2to1_nub, int em2to1_nx, int em2to1_ii, data_type_t *em2to1_out_ub);

void abXYc_matrix2D(data_type_t *abXYcm2D_L, data_type_t *abXYcm2D_Ra, data_type_t *abXYcm2D_Rb, data_type_t abXYcm2D_a, data_type_t abXYcm2D_b, int abXYcm2D_nub, int abXYcm2D_nx, data_type_t abXYcm2D_c);

void abXYc_matrix1D(data_type_t *abXYcm1D_L, data_type_t *abXYcm1D_Ra, data_type_t *abXYcm1D_Rb, data_type_t abXYcm1D_a, data_type_t abXYcm1D_b, int abXYcm1D_nx, data_type_t abXYcm1D_c);

void cshift_matrix2D_dim2(data_type_t *cm2d2_ub, int cm2d2_nub, int cm2d2_nx, int cm2d2_logical, data_type_t *cm2d2_out_ub);

void limiter_matrix2D(data_type_t *lm2D_ubA, data_type_t *lm2D_ubB, int lm2D_nub, int lm2D_nx, data_type_t *lm2D_out_ub);

#endif
