#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#include <stdio.h>
#include <stdlib.h>

#include "../data_type.h"
#include "../array_definition.h"

#include "spu_function.h"
#include "value_spu.h"
#include "barrier.h"

// Local store structures and buffers.
volatile variable_context ctx;
volatile variable_context ctx_update;

int main(unsigned long long spu_id __attribute__ ((unused)), unsigned long long main_data, unsigned long long update_data)
{
//	prepare for the barrier function
//      get index first
  int SPE_id,out_sig;
  SPE_id = spu_read_in_mbox();

 // Send message to PPE asking for notify area address
  spu_write_out_mbox(1);
 //set an unique signal value depending on the SPE ID, we basically assign a different bit to every SPE
   signal[3]=1<<SPE_id;
 //compute the sum of every SPE signal (it will be used to check if every SPE has sent his message)
   total=(1<<MAX_SPE_THREADS)-1;
 //create the array of adresses
   addr=(unsigned int*)malloc(MAX_SPE_THREADS*2*sizeof(unsigned int));
  // Wait for PPE to send the address of the array
  while (spu_stat_in_mbox () < 1);
  //get the address
  tabAddr = spu_read_in_mbox();
  //DMA the array
  mfc_get(addr, tabAddr, MAX_SPE_THREADS*2*sizeof(unsigned int), 31, 0, 0);
  mfc_write_tag_mask(1<<31);
  mfc_read_tag_status_all();

//	finished preparing

  //volatile data_type_t tmp1;
  volatile data_type_t tmp1[2] ;

  unsigned int tag_id;

//	Reserve a tag ID
  tag_id = mfc_tag_reserve();

  spu_writech(MFC_WrTagMask, -1);
//	Input parameter data is a pointer to the data context.
//	Fetch the context, waiting for it to complete.
//	main data matrix
  spu_mfcdma32((void *)(&ctx), (unsigned int)main_data, sizeof(variable_context), tag_id, MFC_GET_CMD);
  spu_mfcdma32((void *)(&ctx_update), (unsigned int)update_data, sizeof(variable_context), tag_id, MFC_GET_CMD);
  (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);
//
  int iter; 
  data_type_t t, tf; 

  iter=0;
  t=0E0;
  tf=ctx.nx*4E0;

do {
//      just output to check result
if ((ctx.Y_location==0)&&(ctx.Z_location==0))
{
                spu_mfcdma32((void *)(&tmp1),
                (unsigned int)(ctx.u+matrix4D(5,(box_nx),(box_ny),(box_nz),(5-1),(1-1),(1-1),(1-1))),
                2*sizeof(data_type_t), tag_id, MFC_GETB_CMD);
                (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);
//
if (sizeof(data_type_t)==4)
{
printf("t=      %E,     %i,     %E,     %E\n",t,iter,tmp1[0],tmp1[1]);
}
else if (sizeof(data_type_t)==8)
{
printf("t=      %16.13E,     %i,     %17.15E,     %17.15E\n",t,iter,tmp1[0],tmp1[1]);
}
}

  iter=iter+1;
  //ctx.dt=0.25E0
  ctx.dt=0.9E0*calcfl(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.Y_location,&ctx.Z_location,SPE_id,tag_id,addr);
  if (ctx.dt>(tf-t)/2.0E0) ctx.dt=(tf-t)/2.0E0;
  t=t+2.0E0*ctx.dt;

//	done output
//----------------------------------------
//----------------------------------------
//	NOTE!
//	THIS VERSION IS ONLY FOR EQUAL DEMENSION
//	ONLY X=Y=Z CAN WORK FOR THIS CODE
//----------------------------------------
//----------------------------------------
  barrier(SPE_id);
  fluidx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,&ctx.Y_location,&ctx.Z_location,tag_id);
  barrier(SPE_id);
  advectbyzx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,ctx.adv_tmp,&ctx.Y_location,&ctx.Z_location,tag_id,SPE_id);
  barrier(SPE_id);
//	the y sweep
  transpose12(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b,&ctx.Y_location,&ctx.Z_location,tag_id,SPE_id);
  barrier(SPE_id);
  copy_matrix(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b,&ctx.Y_location,&ctx.Z_location,tag_id);
  barrier(SPE_id);
  fluidx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,&ctx.Y_location,&ctx.Z_location,tag_id);
  barrier(SPE_id);
  advectbyzx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,ctx.adv_tmp,&ctx.Y_location,&ctx.Z_location,tag_id,SPE_id);
  barrier(SPE_id);

//	z sweep
  transpose13(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b,&ctx.Y_location,&ctx.Z_location,tag_id,SPE_id);
  barrier(SPE_id);
  copy_matrix(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b,&ctx.Y_location,&ctx.Z_location,tag_id);
  barrier(SPE_id);
  fluidx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,&ctx.Y_location,&ctx.Z_location,tag_id);
  barrier(SPE_id);
  advectbyzx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,ctx.adv_tmp,&ctx.Y_location,&ctx.Z_location,tag_id,SPE_id);
  barrier(SPE_id);
  advectbyzx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,ctx.adv_tmp,&ctx.Y_location,&ctx.Z_location,tag_id,SPE_id);
  barrier(SPE_id);
  fluidx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,&ctx.Y_location,&ctx.Z_location,tag_id);
  barrier(SPE_id);

//	back

  transpose13(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b,&ctx.Y_location,&ctx.Z_location,tag_id,SPE_id);
  barrier(SPE_id);
  copy_matrix(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b,&ctx.Y_location,&ctx.Z_location,tag_id);
  barrier(SPE_id);
  advectbyzx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,ctx.adv_tmp,&ctx.Y_location,&ctx.Z_location,tag_id,SPE_id);
  barrier(SPE_id);
  fluidx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,&ctx.Y_location,&ctx.Z_location,tag_id);
  barrier(SPE_id);

//	x again
  transpose12(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b,&ctx.Y_location,&ctx.Z_location,tag_id,SPE_id);
  barrier(SPE_id);
  copy_matrix(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,ctx_update.u,ctx_update.b,&ctx.Y_location,&ctx.Z_location,tag_id);
  barrier(SPE_id);
  advectbyzx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,ctx.adv_tmp,&ctx.Y_location,&ctx.Z_location,tag_id,SPE_id);
  barrier(SPE_id);
  fluidx(ctx.u,ctx.b,&ctx.nx,&ctx.ny,&ctx.nz,&ctx.dt,&ctx.Y_location,&ctx.Z_location,tag_id);
  barrier(SPE_id);
//
   } while (t<tf);

  return (0);

}
