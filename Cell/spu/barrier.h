#ifndef BARRIER_H 
#define BARRIER_H 

//tag mask associated with signals
#define mask_signal 19

// Local store structures and buffers.
volatile variable_context ctx;
volatile variable_context ctx_update;

unsigned int signal[4] __attribute__ ((aligned(128)));
char * ls=((char*)&signal[0])+12;
unsigned int *addr;
//
unsigned int check;
//address in main memory of the array of addresses
unsigned int tabAddr;

unsigned int total;

/*
 * Send the SPE signal to addr
 * addr must be a valid signal notify address
 */
inline void send_signal(unsigned int addr){

        mfc_sndsig(ls,addr,mask_signal,0,0);

}

/*
 * Wait for a signal and read it
 * we don't really care about the value of the signal
 */

inline int receive_signal(){
        while(!spu_stat_signal1());
        return spu_read_signal1();
}

/*
 * The barrier function
 * it must be called by every SPEs
 * SPEs send a signal to SPE0 to say they are ready and then wait for an answer
 * SPE0 waits every signal to be sent and then starts sending signal to unblock the others
 * we use a classic hypercube communication scheme to wake up every SPE
 */

void barrier(int SPE_id)
{

  int i;
  int temp=1;

  if(SPE_id!=0){
        send_signal(addr[0]);
        receive_signal();
  }

  else{
        check = 1;
       while(check < total){
           /*between each call to spu_read_signal1 the notify area is reset to 0
            *since we have used the SPE_CFG_SIGNOTIFY1_OR flag, we are sure that the notify area has not be
            *overwritten by multiple SPEs between 2 calls but has been bit wise OR.
            *each SPE is assigned a different bit so we know we are done when check is equal 11, 1111, 11111111
            *or 1111111111111111 depending on the dimension of the hypercube.
            */

           check|=spu_read_signal1();
       }
  }
  //wake up every SPE using a classic hypercube communication scheme

  for(i=0;i<hypercubeD;i++){

        if(SPE_id<temp){
                send_signal(addr[SPE_id+temp]);
        }
        temp=temp<<1;
  }
}

#endif
