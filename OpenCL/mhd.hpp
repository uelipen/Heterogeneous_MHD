#ifndef MHD_DECLARATION_H_
#define MHD_DECLARATION_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <CL/cl.h>

#include "parameter.hpp"
/*
//
    extern cl_device_type mhdDeviceType;	// device type
    extern cl_context mhdComputeContext;      	// compute context
    extern cl_device_id *mhdComputeDevice;     	// compute device id 
    extern cl_command_queue mhdCommandQueue;   	// compute command queue
//
    extern cl_program fluidx_Program; 	        // compute program
    extern cl_kernel fluidx_Kernel;   		// compute kernel
//
    extern size_t mhdMaxWorkGroupSize;       	// Max allowed work-items in a group
    extern cl_uint mhdMaxDimensions;          	// Max group dimensions allowed 
    extern size_t *mhdMaxWorkItemSizes;      	// Max work-items sizes in each dimensions
    extern cl_ulong mhdTotalLocalMemory;      	// Max local memory allowed
    extern cl_ulong mhdUsedLocalMemory;       	// Used local memory
    extern size_t mhdKernelWorkGroupSize;	// Kernel work group size
//
    extern cl_mem u;                		// fluid component in device memory
    extern cl_mem b;                		// magnetic component in device memory
    extern cl_int nx;				// x length in device memory
    extern cl_int ny;				// y length in device memory
    extern cl_int nz;				// z length in device memory
    extern cl_float dt;				// cfl time in device memory
    extern cl_mem adv_tmpB;			// temporary matrix for b update in device memory
    extern cl_mem u_update;			// fluid update component in device memory
    extern cl_mem b_update;			// magnetic update component in device memory
//
    extern data_type_t *init_u;			// initialization of u in host memory
    extern data_type_t *init_b;			// initialization of b in host memory
    extern int *init_nx;			// initialization of nx in host memory
    extern int *init_ny;			// initialization of ny in host memory
    extern int *init_nz;			// initialization of nz in host memory
    extern data_type_t *init_dt;		// initialization of dt in host memory
    extern data_type_t *init_adv_tmpB;		// initialization of adv_tmpB in host memory
    extern data_type_t *init_u_update;		// initialization of u_update component in host memory
    extern data_type_t *init_b_update;		// initialization of b_update component in host memory
//
*/
/**
 * MHD 
 * Class implements OpenCL MHD sample
 */

class MHD
{

private:

public:

//    ~MHD();

    /**
 *      * initialize the value of matrix u and b
 *           *
 *                */
    void init_value(data_type_t *initV_u, data_type_t *initV_b, int *initV_nx, int *initV_ny, int *initV_nz, data_type_t *init_adv_tmpB, data_type_t *init_u_update, data_type_t *init_b_update);
    /**
     * Allocate and initialize value in host memory
     */
    int setupMHD_initialization();

    /**
     * OpenCL related initialisations. 
     * Set up Context, Device list, Command Queue, Memory buffers
     */
    int setupCL_ContextCommandMemory();

    /**
     * Build CL Program and Kernel executable 
     */
    int setupCL_ProgramKernel();

    /**
     * Set values for kernels' arguments
     */
    int setupCL_ArgumentKernel();

    /**
     * Enqueue calls to the kernels
     * on to the command queue, wait till end of kernel execution.
     */
    int runCL_Kernel();

    /**
     * Combine all the setup sub function above
     * of execution domain, perform all sample setup
     */
    int setup();

    /**
     * Combine all the run sub function above 
     * Run OpenCL MHD
     */
    int run();

    /**
     * Cleanup memory allocations
     */
    int cleanup();

};

#endif // MHD_DECLARATION_H_ 
