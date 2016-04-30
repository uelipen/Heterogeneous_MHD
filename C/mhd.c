//	TVD split MHD code on CELL BE
//	copyright (C) 2009, Bijia PANG
//	written August 2009 by Bijia PANG, bpang@physics.utoronto.ca
//	
//
#include <stdio.h>
#include <stdlib.h>

#include <libspe2.h>
#include <pthread.h>

#include "parameter.h"
#include "array_definition.h"
#include "ppu_function.h"
#include "spu_function.h"
#include "value_spu.h"

//#include "init.c"

float main_u[nu*TOTAL_CELL_NUMBER];// __attribute__ ((aligned (16)));
float main_b[nb*TOTAL_CELL_NUMBER];// __attribute__ ((aligned (16)));

float update_u[nu*TOTAL_CELL_NUMBER];// __attribute__ ((aligned (16)));
float update_b[nb*TOTAL_CELL_NUMBER];// __attribute__ ((aligned (16)));

typedef struct ppu_pthread_data {
  spe_context_ptr_t spe_mhd;
  pthread_t pthread;
  void *argp;
  void *envp;
} ppu_pthread_data_t;

int main()
{
printf("program start!\n");

  ppu_pthread_data_t data;
//	main data box
  variable_context ctx;
//	update data box
  variable_context ctx_update;

  /* Initialize context run data */
//      call init
  ctx.u = main_u;
  ctx.b = main_b;
  ctx_update.u = update_u;
  ctx_update.b = update_b;

  init(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.total_cell_number,&ctx.no_meaning);
  init(ctx_update.u,ctx_update.b,&ctx_update.nx,&ctx_update.ny,&ctx_update.nz,&ctx_update.total_cell_number,&ctx_update.no_meaning);

int i,j,k,ii;
  float t, tf;
  int iter;
  iter=0;
  t=0;
  tf=ctx.nx*4;
do {
  iter=iter+1;
  //ctx.dt=0.5;//0.9*cfl(u,b)
  ctx.dt=0.9*calcfl(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz);
  if (ctx.dt>(tf-t)/2) ctx.dt=(tf-t)/2;
  t=t+2*ctx.dt;

printf("t=      %e,     %i,     %e,     %e\n",t,iter,ctx.u[matrix4D(5,ctx.nx,ctx.ny,ctx.nz,(5-1),(1-1),(3-1),(1-1))],ctx.u[matrix4D(5,ctx.nx,ctx.ny,ctx.nz,(1-1),(2-1),(1-1),(1-1))]);


//	printf("----------- %e\n",ctx.dt);
  //check_value(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b,iter);
  fluidx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,ctx_update.u,ctx_update.b,iter);
  copy_matrix(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b);
  //check_update(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b,iter);
  advectbyzx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,ctx_update.u,ctx_update.b);
  copy_matrix(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b);

//	the y sweep
  transpose12(ctx.u,ctx.b,ctx_update.u,ctx_update.b,&ctx.nx,&ctx.ny,&ctx.nz,iter);
  copy_matrix(ctx.u,ctx.b,&ctx.ny,&ctx.nx,&ctx.nz,ctx_update.u,ctx_update.b);
  fluidx(ctx.u,ctx.b,&ctx.ny,&ctx.nx,&ctx.nz,&ctx.dt,ctx_update.u,ctx_update.b,iter);
  copy_matrix(ctx.u,ctx.b,&ctx.ny,&ctx.nx,&ctx.nz,ctx_update.u,ctx_update.b);
  check_value(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b,iter);
  //check_update(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b,iter);
  advectbyzx(ctx.u,ctx.b,&ctx.ny,&ctx.nx,&ctx.nz,&ctx.dt,ctx_update.u,ctx_update.b);
  copy_matrix(ctx.u,ctx.b,&ctx.ny,&ctx.nx,&ctx.nz,ctx_update.u,ctx_update.b);
//	z sweep

  transpose13(ctx.u,ctx.b,ctx_update.u,ctx_update.b,&ctx.ny,&ctx.nx,&ctx.nz);
  copy_matrix(ctx.u,ctx.b,&ctx.nz,&ctx.nx,&ctx.ny,ctx_update.u,ctx_update.b);

  fluidx(ctx.u,ctx.b,&ctx.nz,&ctx.nx,&ctx.ny,&ctx.dt,ctx_update.u,ctx_update.b,iter);
  copy_matrix(ctx.u,ctx.b,&ctx.nz,&ctx.nx,&ctx.ny,ctx_update.u,ctx_update.b);
  advectbyzx(ctx.u,ctx.b,&ctx.nz,&ctx.nx,&ctx.ny,&ctx.dt,ctx_update.u,ctx_update.b);
  copy_matrix(ctx.u,ctx.b,&ctx.nz,&ctx.nx,&ctx.ny,ctx_update.u,ctx_update.b);
  advectbyzx(ctx.u,ctx.b,&ctx.nz,&ctx.nx,&ctx.ny,&ctx.dt,ctx_update.u,ctx_update.b);
  copy_matrix(ctx.u,ctx.b,&ctx.nz,&ctx.nx,&ctx.ny,ctx_update.u,ctx_update.b);
  fluidx(ctx.u,ctx.b,&ctx.nz,&ctx.nx,&ctx.ny,&ctx.dt,ctx_update.u,ctx_update.b,iter);
  copy_matrix(ctx.u,ctx.b,&ctx.nz,&ctx.nx,&ctx.ny,ctx_update.u,ctx_update.b);

//	back

  transpose13(ctx.u,ctx.b,ctx_update.u,ctx_update.b,&ctx.nz,&ctx.nx,&ctx.ny);
  copy_matrix(ctx.u,ctx.b,&ctx.ny,&ctx.nx,&ctx.nz,ctx_update.u,ctx_update.b);

  advectbyzx(ctx.u,ctx.b,&ctx.ny,&ctx.nx,&ctx.nz,&ctx.dt,ctx_update.u,ctx_update.b);
  copy_matrix(ctx.u,ctx.b,&ctx.ny,&ctx.nx,&ctx.nz,ctx_update.u,ctx_update.b);
  fluidx(ctx.u,ctx.b,&ctx.ny,&ctx.nx,&ctx.nz,&ctx.dt,ctx_update.u,ctx_update.b,iter);
  copy_matrix(ctx.u,ctx.b,&ctx.ny,&ctx.nx,&ctx.nz,ctx_update.u,ctx_update.b);

//	x again

  transpose12(ctx.u,ctx.b,ctx_update.u,ctx_update.b,&ctx.ny,&ctx.nx,&ctx.nz,iter);
  copy_matrix(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b);

  advectbyzx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,ctx_update.u,ctx_update.b);
  copy_matrix(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b);
  fluidx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,ctx_update.u,ctx_update.b,iter);
  copy_matrix(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b);


//printf("t=	%e,	%i,	%e,	%e\n",t,iter,ctx.u[matrix4D(5,ctx.nx,ctx.ny,ctx.nz,(5-1),(1-1),(1-1),(1-1))],ctx.u[matrix4D(5,ctx.nx,ctx.ny,ctx.nz,(1-1),(1-1),(1-1),(1-1))]);
	
   //} while (iter<1);
   } while (t<tf);

printf("program finished!\n");
  return (0);

}


