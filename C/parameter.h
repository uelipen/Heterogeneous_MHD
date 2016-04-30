//	parameters for simulation

#ifndef PARAMETER_H
#define PARAMETER_H

#define n 100
#define box_nx n
#define box_ny n
#define box_nz n
#define TOTAL_CELL_NUMBER box_nx*box_ny*box_nz
#define nu 5
#define nb 3

typedef struct {
  float *u;
  float *b;
  int nx;
  int ny;
  int nz;
  float dt;
  int total_cell_number;
  int no_meaning;
} variable_context;

#endif
