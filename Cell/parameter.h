//	parameters for simulation

#ifndef PARAMETER_H
#define PARAMETER_H

#include "data_type.h"

#define n 16
#define box_nx n
#define box_ny n
#define box_nz n
#define TOTAL_CELL_NUMBER (box_nx*box_ny*box_nz)
#define nu 5
#define nb 3
#define size_X_trans_3D 4
#define size_Y_trans_3D 4
#define size_Z_trans_3D 4
#define size_trans_3D (size_X_trans_3D*size_Y_trans_3D*size_Z_trans_3D)

#define MAX_SPE_THREADS 16
#define num_SPE_Y 4
#define num_SPE_Z 4
#define hypercubeD 4

#endif
