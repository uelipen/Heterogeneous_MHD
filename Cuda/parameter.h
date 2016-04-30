#ifndef PARAMETER_H 
#define PARAMETER_H

#define n 128//64//32//16//128
#define box_nx n
#define box_ny n
#define box_nz n
#define total_cell_number box_nx*box_ny*box_nz
#define BLOCK_SIZE box_nx
#define nu 5
#define nb 3

typedef struct {
  float *u;
  float *b;
  int nx;
  int ny;
  int nz;
  float dt;
} variable_context;

#endif
