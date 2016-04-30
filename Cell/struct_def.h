//	struct definition for simulation

#ifndef STRUCT_DEF_H
#define STRUCT_DEF_H

typedef struct{
  float *u;
  float *b;
  int nx;
  int ny;
  int nz;
  float *adv_tmp;
  float *trans_u_3D;
  float *trans_b_3D;
//
  float dt;
  int Y_location;
  int Z_location;
  float *reserved_1;
}variable_context_float;

typedef struct{
  double *u;
  double *b;
  int nx;
  int ny;
  int nz;
  double *adv_tmp;
  double *trans_u_3D;
  double *trans_b_3D;
//
  double dt;
  int Y_location;
  int Z_location;
}variable_context_double;

#endif
