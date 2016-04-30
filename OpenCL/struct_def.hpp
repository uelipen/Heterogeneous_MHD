#ifndef STRUCT_DEF_H_
#define STRUCT_DEF_H_

#include "data_type.hpp"

typedef struct{
  data_type_t *u;		// fluid
  data_type_t *b;		// magnetic
  int nx;			// X length
  int ny;			// Y length
  int nz;			// Z length
  data_type_t dt;		// cfl time
  data_type_t *adv_tmpB;	// matrix for magnetic in magnetic update function
  data_type_t *u_update;	// matrix for fluid in transpose function
  data_type_t *b_update;	// matrix for magnetic in transpose function
}variable_context_host;

typedef struct{
  cl_mem u;                  	// fluid
  cl_mem b;                     // magnetic
  cl_int nx;                    // X length
  cl_int ny;                    // Y length
  cl_int nz;                    // Z length
  cl_float dt;                  // cfl time
  cl_mem adv_tmpB;              // matrix for magnetic in magnetic update function
  cl_mem u_update;              // matrix for fluid in transpose function
  cl_mem b_update;              // matrix for magnetic in transpose function 
}variable_context_device;

#endif	// STRUCT_DEF_H_

