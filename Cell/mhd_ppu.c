//	TVD split MHD code on CELL BE
//	copyright (C) 2009, Bijia PANG
//	written August 2009 by Bijia PANG, bpang@physics.utoronto.ca
//	
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <libspe2.h>
#include <pthread.h>

#include "data_type.h"

#include "parameter.h"
#include "array_definition.h"
#include "ppu_function.h"

void PPE_check_value(data_type_t *PPE_checkV_u, data_type_t *PPE_checkV_b, int *PPE_checkV_nx, int *PPE_checkV_ny, int *PPE_checkV_nz);

data_type_t main_u[nu*TOTAL_CELL_NUMBER] __attribute__ ((aligned (16)));
data_type_t main_b[nb*TOTAL_CELL_NUMBER] __attribute__ ((aligned (16)));

data_type_t main_u_update[nu*TOTAL_CELL_NUMBER] __attribute__ ((aligned (16)));
data_type_t main_b_update[nb*TOTAL_CELL_NUMBER] __attribute__ ((aligned (16)));

data_type_t main_adv_tmp[TOTAL_CELL_NUMBER] __attribute__ ((aligned (16)));

data_type_t transpose_u_3D[5*size_trans_3D] __attribute__ ((aligned (16)));
data_type_t transpose_b_3D[3*size_trans_3D] __attribute__ ((aligned (16)));

extern spe_program_handle_t mhd_spu;

typedef struct ppu_pthread_data {
  spe_context_ptr_t spe_mhd;
  pthread_t pthread;
  void *argp;
  void *envp;
} ppu_pthread_data_t;

void *ppu_pthread_function(void *arg) {
  ppu_pthread_data_t *datap = (ppu_pthread_data_t *)arg;
  unsigned int entry = SPE_DEFAULT_ENTRY;
  if (spe_context_run(datap->spe_mhd, &entry, 0, datap->argp, datap->envp, NULL) < 0) {
    perror ("Failed running context");
    exit (1);
  }
  pthread_exit(NULL);
}

int main()
{
printf("program start!\n");

  ppu_pthread_data_t data_all __attribute__ ((aligned (16))); 
  ppu_pthread_data_t datas[MAX_SPE_THREADS] __attribute__ ((aligned (16))); 
//	main data box
  variable_context ctx_all __attribute__ ((aligned (16)));
  variable_context ctxs[MAX_SPE_THREADS] __attribute__ ((aligned (16)));
  variable_context ctx_update_all __attribute__ ((aligned (16)));
  variable_context ctx_updates[MAX_SPE_THREADS] __attribute__ ((aligned (16)));

//	Determine the number of SPE threads to create.
  int spe_threads;
  spe_threads = spe_cpu_info_get(SPE_COUNT_USABLE_SPES, -1);
  if (spe_threads > MAX_SPE_THREADS) spe_threads = MAX_SPE_THREADS;

  /* Initialize context run data */
//      call init to initilize the whole box
  ctx_all.u = main_u;
  ctx_all.b = main_b;
  ctx_update_all.u = main_u_update;
  ctx_update_all.b = main_b_update;

  ctx_all.adv_tmp = main_adv_tmp;

  init(ctx_all.u,ctx_all.b,&ctx_all.nx,&ctx_all.ny,&ctx_all.nz,ctx_all.adv_tmp);

  ctx_all.trans_u_3D = transpose_u_3D;
  ctx_all.trans_b_3D = transpose_b_3D;
//
//	set value for each SPE
//	4 SPE, 2*2 for Y*Z plane
//	16 SPE, 4*4 for Y*Z plane
  int Y_inter=box_ny/num_SPE_Y;
  int Z_inter=box_nz/num_SPE_Z;
  int index_Y,index_Z;
  int index_spe_u,index_spe_b,index_spe_adv_tmp;
  int i;
  if (((box_ny%num_SPE_Y)!=0)||((box_nz%num_SPE_Z)!=0)) 
  {
	printf("error: box size can't be divided by SPE number!\n");
	exit (1);
  }
  if (MAX_SPE_THREADS!=(num_SPE_Y*num_SPE_Z))
  {
	printf("error: MAX_SPE_THREADS is NOT equal to num_SPE_Y*num_SPE_Z");
	exit (1);
  }
  if (num_SPE_Y!=num_SPE_Z)
  {
	printf("error: num_SPE_Y must equal to num_SPE_Z for transpose matrix");
	exit (1);
  }
  for (i=0; i<spe_threads; i++) 
  {
	index_Y=(i%num_SPE_Y);
	index_Z=(i/num_SPE_Y);
	index_spe_u=matrix4D(5,box_nx,box_ny,box_nz,(1-1),(1-1),(index_Y*Y_inter),(index_Z*Z_inter));
	index_spe_b=matrix4D(3,box_nx,box_ny,box_nz,(1-1),(1-1),(index_Y*Y_inter),(index_Z*Z_inter));
	index_spe_adv_tmp=matrix3D(box_nx,box_ny,box_nz,(1-1),(index_Y*Y_inter),(index_Z*Z_inter));
//
	ctxs[i].u = ctx_all.u;
	ctxs[i].b = ctx_all.b;
	ctxs[i].nx = ctx_all.nx;
	ctxs[i].ny = ctx_all.ny/num_SPE_Y;
	ctxs[i].nz = ctx_all.nz/num_SPE_Z;
	ctx_updates[i].u = ctx_update_all.u;
	ctx_updates[i].b = ctx_update_all.b;
	ctxs[i].adv_tmp = ctx_all.adv_tmp;
	ctxs[i].Y_location = index_Y; 
	ctxs[i].Z_location = index_Z;
  }

  for (i=0; i<spe_threads; i++)
  {
  /* Create a SPE context */
  //if ((datas[i].spe_mhd = spe_context_create (0, NULL)) == NULL) {
 	if (i==0)
	{
  if ((datas[0].spe_mhd = spe_context_create (SPE_CFG_SIGNOTIFY1_OR | SPE_MAP_PS, NULL)) == NULL) {
    perror ("Failed creating context");
    exit (1);
  }
	}
	else
	{
  if ((datas[i].spe_mhd = spe_context_create (SPE_MAP_PS, NULL)) == NULL) {
    perror ("Failed creating context");
    exit (1);
  }
	}
  /* Load SPE program into the SPE context*/
  if (spe_program_load (datas[i].spe_mhd, &mhd_spu))  {
    perror ("Failed loading program");
    exit (1);
  }
//
  datas[i].argp=&ctxs[i];
  datas[i].envp=&ctx_updates[i];
//
  /* Create pthread for each of the SPE contexts */
  if (pthread_create (&datas[i].pthread, NULL, &ppu_pthread_function, &datas[i])) {
    perror ("Failed creating thread");
    exit (1);
  }
  }
//	send SPE id to each SPE using mailbox
  for (i=0; i<spe_threads; i++)
  {
        spe_in_mbox_write(datas[i].spe_mhd,(unsigned int*)&i,1,SPE_MBOX_ANY_NONBLOCKING);
  }
        printf("PPE sending adresses\n");

                //pointer to the signal notify area 1
                spe_sig_notify_1_area_t *test;
                //pointer to the LS area (not useful to call barrier but you might need it if you want to DMA between SPEs)
                void * test2;
                //array of adresses (notify area and LS adresses of each SPE)
                unsigned int *addr;

                //we align it on a 128 bits basis since it will be DMAed to each SPE
                addr=(unsigned int*)malloc(127+2*MAX_SPE_THREADS*sizeof(unsigned int));
                while ((int)(addr) & 0x7f) addr++;

                //we can store the message receive from each SPE in this variable
                unsigned int message[MAX_SPE_THREADS];

         //get addresses
         for (i=0; i<MAX_SPE_THREADS; ++i){
                 //we get a pointer to a struct containing the data we want
                 test=(spe_sig_notify_1_area_t *)spe_ps_area_get(datas[i].spe_mhd,SPE_SIG_NOTIFY_1_AREA);
                 //we get the real address of the notify area 1
                 addr[i]=(unsigned int) &(test->SPU_Sig_Notify_1);
                 //we get the LS area address
                 test2=spe_ls_area_get(datas[i].spe_mhd);
	 	addr[MAX_SPE_THREADS+i]=(unsigned int)test2;
         }
         // We wait for each SPE 's mail and we send them the address of the array.

         for (i=0; i<MAX_SPE_THREADS; ++i){
             while (spe_out_mbox_status(datas[i].spe_mhd) <= 0);
             spe_out_mbox_read(datas[i].spe_mhd, &message[i], 1);
             spe_in_mbox_write(datas[i].spe_mhd, (void*)&addr, 1, SPE_MBOX_ANY_NONBLOCKING);
         }
printf("finished signal stuff\n");
//	finished signal stuff
//
  for (i=0; i<spe_threads; i++)
  {
  /* Wait for the threads to complete */
  if (pthread_join (datas[i].pthread, NULL)) {
    perror ("Failed joining thread\n");
    exit (1);
  }
  /* Destroy SPE context */
  if (spe_context_destroy (datas[i].spe_mhd) != 0) {
    perror("Failed destroying context");
    exit (1);
  }
  }

//PPE_check_value(ctx_all.u, ctx_all.b, &ctx_all.nx, &ctx_all.ny, &ctx_all.nz);

printf("program finished!\n");
  return (0);

}


