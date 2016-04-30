#ifndef VALUE_SPU_H
#define VALUE_SPU_H 

#include "parameter.h"

#define u_memSize_WHOLE 5*box_nx*box_ny*box_nz
#define b_memSize_WHOLE 3*box_nx*box_ny*box_nz
#define CELL_PER_BLOCK box_nx
#define u_memSize_PER_BLOCK 5*CELL_PER_BLOCK
#define b_memSize_PER_BLOCK 3*CELL_PER_BLOCK

#endif

