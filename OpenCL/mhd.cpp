#include <fstream>
#include <fcntl.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <stdlib.h>
#include <malloc.h>
#include <SDKUtil/SDKFile.hpp>
#include <CL/cl.h>

#include "data_type.hpp"
#include "struct_def.hpp"
#include "info_def.hpp"
#include "parameter.hpp"
#include "array_definition.hpp"
#include "mhd.hpp"

//
    cl_device_type mhdDeviceType;        // device type
    cl_context mhdComputeContext;        // compute context
    cl_device_id *mhdComputeDevice;      // compute device id
    cl_command_queue mhdCommandQueue;    // compute command queue
    cl_program fluidx_Program;           // compute program
    cl_kernel fluidx_Kernel;             // compute kernel
    size_t mhdMaxWorkGroupSize;          // Max allowed work-items in a group
    cl_uint mhdMaxDimensions;            // Max group dimensions allowed
    size_t *mhdMaxWorkItemSizes;         // Max work-items sizes in each dimensions
    cl_ulong mhdTotalLocalMemory;        // Max local memory allowed
    cl_ulong mhdUsedLocalMemory;         // Used local memory
    size_t mhdKernelWorkGroupSize;       // Kernel work group size

    variable_context_host ctx_host;
    variable_context_device ctx_device;

//
int
MHD::setupMHD_initialization()
{

    //init_u = (data_type_t*)malloc(total_cell_number*nu*sizeof(data_type_t));
    ctx_host.u = (data_type_t*)memalign(16,total_cell_number*nu*sizeof(data_type_t));
    if(ctx_host.u == NULL)
    {
	printf("Error: Failed to allocate ctx_host.u host memory!\n");    
        return MHD_FAILURE;
    }
//
    ctx_host.b = (data_type_t*)memalign(16,total_cell_number*nb*sizeof(data_type_t));
    if(ctx_host.b == NULL)
    {
	printf("Error: Failed to allocate ctx_host.b host memory!\n");
        return MHD_FAILURE;
    }
//
    ctx_host.adv_tmpB = (data_type_t*)memalign(16,total_cell_number*sizeof(data_type_t));
    if(ctx_host.adv_tmpB == NULL)
    {
        printf("Error: Failed to allocate ctx_host.adv_tmpB host memory!\n");
        return MHD_FAILURE;
    }
//
    ctx_host.u_update = (data_type_t*)memalign(16,total_cell_number*nu*sizeof(data_type_t));
    if(ctx_host.u_update == NULL)
    {
        printf("Error: Failed to allocate ctx_host.u_update host memory!\n");
        return MHD_FAILURE;
    }
//
    ctx_host.b_update = (data_type_t*)memalign(16,total_cell_number*nb*sizeof(data_type_t));
    if(ctx_host.b_update == NULL)
    {
        printf("Error: Failed to allocate ctx_host.b_update host memory!\n");
        return MHD_FAILURE;
    }
//
    MHD::init_value(ctx_host.u,ctx_host.b,&ctx_host.nx,&ctx_host.ny,&ctx_host.nz,ctx_host.adv_tmpB,ctx_host.u_update,ctx_host.b_update);
// 
	return MHD_SUCCESS;
}

void 
MHD::init_value(data_type_t *initV_u, data_type_t *initV_b, int *initV_nx, int *initV_ny, int *initV_nz, data_type_t *initV_adv_tmpB, data_type_t *initV_u_update, data_type_t *initV_b_update)
{
printf("beginning in init\n");

  *initV_nx = box_nx;
  *initV_ny = box_ny;
  *initV_nz = box_nz;

  int i,j,k,ii;
  FILE *initV_File;
  initV_File=fopen("source_init_alfvenlinear_16cube.dat","r+");

  for (i=0;i<(*initV_nx);i++)
  {
        for (j=0;j<(*initV_ny);j++)
        {
                for (k=0;k<(*initV_nz);k++)
                {
                        initV_adv_tmpB[a3D_FinC((*initV_nx),(*initV_ny),(*initV_nz),i,j,k)]=0.0E0;
for (ii=0;ii<5;ii++)
{
                        initV_u_update[a4D_FinC(5,(*initV_nx),(*initV_ny),(*initV_nz),ii,i,j,k)]=0.0E0;
}
for (ii=0;ii<3;ii++)
{
                        initV_b_update[a4D_FinC(3,(*initV_nx),(*initV_ny),(*initV_nz),ii,i,j,k)]=0.0E0;
}

                        if (sizeof(data_type_t)==4)
                        {
                        fscanf(initV_File,"%f",&initV_u[a4D_FinC(5,(*initV_nx),(*initV_ny),(*initV_nz),0,i,j,k)]);
                        fscanf(initV_File,"%f",&initV_u[a4D_FinC(5,(*initV_nx),(*initV_ny),(*initV_nz),1,i,j,k)]);
                        fscanf(initV_File,"%f",&initV_u[a4D_FinC(5,(*initV_nx),(*initV_ny),(*initV_nz),2,i,j,k)]);
                        fscanf(initV_File,"%f",&initV_u[a4D_FinC(5,(*initV_nx),(*initV_ny),(*initV_nz),3,i,j,k)]);
                        fscanf(initV_File,"%f",&initV_u[a4D_FinC(5,(*initV_nx),(*initV_ny),(*initV_nz),4,i,j,k)]);
                        fscanf(initV_File,"%f",&initV_b[a4D_FinC(3,(*initV_nx),(*initV_ny),(*initV_nz),0,i,j,k)]);
                        fscanf(initV_File,"%f",&initV_b[a4D_FinC(3,(*initV_nx),(*initV_ny),(*initV_nz),1,i,j,k)]);
                        fscanf(initV_File,"%f",&initV_b[a4D_FinC(3,(*initV_nx),(*initV_ny),(*initV_nz),2,i,j,k)]);
                        }
                        else if (sizeof(data_type_t)==8)
                        {
                        fscanf(initV_File,"%lf",&initV_u[a4D_FinC(5,(*initV_nx),(*initV_ny),(*initV_nz),0,i,j,k)]);
                        fscanf(initV_File,"%lf",&initV_u[a4D_FinC(5,(*initV_nx),(*initV_ny),(*initV_nz),1,i,j,k)]);
                        fscanf(initV_File,"%lf",&initV_u[a4D_FinC(5,(*initV_nx),(*initV_ny),(*initV_nz),2,i,j,k)]);
                        fscanf(initV_File,"%lf",&initV_u[a4D_FinC(5,(*initV_nx),(*initV_ny),(*initV_nz),3,i,j,k)]);
                        fscanf(initV_File,"%lf",&initV_u[a4D_FinC(5,(*initV_nx),(*initV_ny),(*initV_nz),4,i,j,k)]);
                        fscanf(initV_File,"%lf",&initV_b[a4D_FinC(3,(*initV_nx),(*initV_ny),(*initV_nz),0,i,j,k)]);
                        fscanf(initV_File,"%lf",&initV_b[a4D_FinC(3,(*initV_nx),(*initV_ny),(*initV_nz),1,i,j,k)]);
                        fscanf(initV_File,"%lf",&initV_b[a4D_FinC(3,(*initV_nx),(*initV_ny),(*initV_nz),2,i,j,k)]);
                        }
                }
        }
  }
  fclose(initV_File);

printf("end in init\n");
}

int
MHD::setupCL_ContextCommandMemory()
{
    cl_int status_in_setupCL = CL_SUCCESS;

    mhdDeviceType = CL_DEVICE_TYPE_GPU;//CL_DEVICE_TYPE_CPU;

    /* Create context from given device type */
    mhdComputeContext = clCreateContextFromType(
		0,
                mhdDeviceType,
                NULL,
                NULL,
                &status_in_setupCL);

    if (status_in_setupCL != CL_SUCCESS)
    {
	printf("Error: Failed to create a compute context!\n");
        return MHD_FAILURE;
    }

    /* First, get the size of device list data */
    size_t deviceListSize;
    status_in_setupCL = clGetContextInfo(
		mhdComputeContext,
                CL_CONTEXT_DEVICES,
                0,
                NULL,
                &deviceListSize);
    if (status_in_setupCL != CL_SUCCESS)
    {
	printf("Error: Failed to get the size of device list!\n");
        return MHD_FAILURE;
    }

    /* Now allocate memory for device list based on the size we got earlier */
    mhdComputeDevice = (cl_device_id *)malloc(deviceListSize);
    if(mhdComputeDevice == NULL) 
    {
	printf("Error: Failed to allocate memory for device list!\n");
        return MHD_FAILURE;
    }

    /* Now, get the device list data */
    status_in_setupCL = clGetContextInfo(
		mhdComputeContext,
                CL_CONTEXT_DEVICES,
                deviceListSize,
                mhdComputeDevice,
                NULL);
    if (status_in_setupCL != CL_SUCCESS)
    {
	printf("Error: Failed to get device list data!\n");
        return MHD_FAILURE;
    }

    /* Create command queue */
    mhdCommandQueue = clCreateCommandQueue(
		mhdComputeContext,
                mhdComputeDevice[0],
                0,
                &status_in_setupCL);
    if (status_in_setupCL != CL_SUCCESS)
    {
	printf("Error: Failed to create command queue!\n");
        return MHD_FAILURE;
    }

    /* Get Device specific Information */
    status_in_setupCL = clGetDeviceInfo(
		mhdComputeDevice[0],
		CL_DEVICE_MAX_WORK_GROUP_SIZE,
                sizeof(size_t),
                (void*)&mhdMaxWorkGroupSize,
                NULL);
printf("Max Work Group Size is %i\n",mhdMaxWorkGroupSize);
    if (status_in_setupCL != CL_SUCCESS)
    {
	printf("Error: Failed to get CL_DEVICE_MAX_WORK_GROUP_SIZE!\n");
        return MHD_FAILURE;
    }

    status_in_setupCL = clGetDeviceInfo(
		mhdComputeDevice[0],
                CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,
                sizeof(cl_uint),
                (void*)&mhdMaxDimensions,
                NULL);
printf("Max Dimensions is %i\n",mhdMaxDimensions);
    if (status_in_setupCL != CL_SUCCESS)
    {
	printf("Error: Failed to get CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS!\n");
        return MHD_FAILURE;
    }

    mhdMaxWorkItemSizes = (size_t*)malloc(mhdMaxDimensions * sizeof(size_t));
    status_in_setupCL = clGetDeviceInfo(
		mhdComputeDevice[0],
                CL_DEVICE_MAX_WORK_ITEM_SIZES,
                sizeof(size_t) * mhdMaxDimensions,
                (void*)mhdMaxWorkItemSizes,
                NULL);
printf("Max Work Item Sizes is %i	%i	%i\n",mhdMaxWorkItemSizes[0],mhdMaxWorkItemSizes[1],mhdMaxWorkItemSizes[2]);
    if (status_in_setupCL != CL_SUCCESS)
    {
	printf("Error: Failed to get CL_DEVICE_MAX_WORK_ITEM_SIZES!\n");
        return MHD_FAILURE;
    }

    status_in_setupCL = clGetDeviceInfo(
		mhdComputeDevice[0],
                CL_DEVICE_LOCAL_MEM_SIZE,
                sizeof(cl_ulong),
                (void *)&mhdTotalLocalMemory,
                NULL);
printf("Total Local Memory is %i\n",mhdTotalLocalMemory);
    if (status_in_setupCL != CL_SUCCESS)
    {
	printf("Error: Failed to get CL_DEVICE_LOCAL_MEM_SIZE!\n");
        return MHD_FAILURE;
    }

    /*
 *      * Create and initialize memory objects
 *           */

    /* Create memory objects */
//	fluid component
    ctx_device.u = clCreateBuffer(
		mhdComputeContext,
		CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
		total_cell_number*nu*sizeof(data_type_t),	
		ctx_host.u,
		&status_in_setupCL);
    if (status_in_setupCL != CL_SUCCESS)
    {
	printf("Error: Failed to allocate fluid component in device memory!\n");
        return MHD_FAILURE;
    }
//	magnetic component
    ctx_device.b = clCreateBuffer(
		mhdComputeContext,
		CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
		total_cell_number*nb*sizeof(data_type_t),
		ctx_host.u,
		&status_in_setupCL);
    if (status_in_setupCL != CL_SUCCESS)
    {
	printf("Error: Failed to allocate magnetic component in device memory!\n");
        return MHD_FAILURE;
    }
//	adv_tmpB
    ctx_device.adv_tmpB = clCreateBuffer(
                mhdComputeContext,
                CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
                total_cell_number*sizeof(data_type_t),
                ctx_host.adv_tmpB,
                &status_in_setupCL);
    if (status_in_setupCL != CL_SUCCESS)
    {
        printf("Error: Failed to allocate adv for temporary in device memory!\n");
        return MHD_FAILURE;
    }
//	u_update
    ctx_device.u_update = clCreateBuffer(
                mhdComputeContext,
                CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
                total_cell_number*nu*sizeof(data_type_t),
                ctx_host.u_update,
                &status_in_setupCL);
    if (status_in_setupCL != CL_SUCCESS)
    {
        printf("Error: Failed to allocate fluid update component in device memory!\n");
        return MHD_FAILURE;
    }
//	b_update
    ctx_device.b_update = clCreateBuffer(
                mhdComputeContext,
                CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
                total_cell_number*nb*sizeof(data_type_t),
                ctx_host.b_update,
                &status_in_setupCL);
    if (status_in_setupCL != CL_SUCCESS)
    {
        printf("Error: Failed to allocate magnetic update component in device memory!\n");
        return MHD_FAILURE;
    }
//
	return MHD_SUCCESS;
}

int
MHD::setupCL_ProgramKernel()
{

    cl_int status_in_setupCL = CL_SUCCESS;

    /* create a CL program using the kernel source */
    streamsdk::SDKFile kernelFile;
    kernelFile.open("fluidx_Kernels.cl");
    const char * source = kernelFile.source().c_str();
    size_t sourceSize[] = { strlen(source) };

    fluidx_Program = clCreateProgramWithSource(
		mhdComputeContext,
	        1,
	        &source,
	        sourceSize,
	        &status_in_setupCL);
    if (status_in_setupCL != CL_SUCCESS)
    {
	printf("Error: Failed to create program for fluidx_Kernels.cl!\n");
        return MHD_FAILURE;
    }

    /* create a cl program executable for all the devices specified */
    status_in_setupCL = clBuildProgram(
		fluidx_Program,
                1,
                &mhdComputeDevice[0],
                NULL,
                NULL,
                NULL);
    if (status_in_setupCL != CL_SUCCESS)
    {
	char buffer[2048];
	size_t len;
	clGetProgramBuildInfo(fluidx_Program, mhdComputeDevice[0], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
	printf("%s\n", buffer);
	printf("Error: Failed to build program executable for fluidx_Kernels.cl!\n");
	return MHD_FAILURE;
    }

    /* get a kernel object handle for a kernel with the given name */
    fluidx_Kernel = clCreateKernel(
		fluidx_Program,
                "fluidx",
                &status_in_setupCL);
    if (status_in_setupCL != CL_SUCCESS)
    {
	printf("Error: Failed to create kernel for fluidx_Kernels.cl!\n");
	return MHD_FAILURE;
    }

        return MHD_SUCCESS;
}

int
MHD::setupCL_ArgumentKernel()
{
    cl_int status_in_setupCL = CL_SUCCESS;

    /*** Set appropriate arguments to the kernel ***/

//	u, fluid component
    status_in_setupCL = clSetKernelArg(
		fluidx_Kernel,
                0,
                sizeof(cl_mem),
                (void *)&ctx_device.u);
    if (status_in_setupCL != CL_SUCCESS)
    {
	printf("Error: Failed to set argument to the kernel! (u)\n");
        return MHD_FAILURE;
    }
//	b, magnetic component
    status_in_setupCL = clSetKernelArg(
                fluidx_Kernel,
                1,
                sizeof(cl_mem),
                (void *)&ctx_device.b);
    if (status_in_setupCL != CL_SUCCESS)
    {
        printf("Error: Failed to set argument to the kernel! (b)\n");
        return MHD_FAILURE;
    }
//	nx, x length
    status_in_setupCL = clSetKernelArg(
                fluidx_Kernel,
                2,
                sizeof(cl_int),
                (void *)&ctx_device.nx);
    if (status_in_setupCL != CL_SUCCESS)
    {
        printf("Error: Failed to set argument to the kernel! (nx)\n");
        return MHD_FAILURE;
    }
//	ny, y length
    status_in_setupCL = clSetKernelArg(
                fluidx_Kernel,
                3,
                sizeof(cl_int),
                (void *)&ctx_device.ny);
    if (status_in_setupCL != CL_SUCCESS)
    {
        printf("Error: Failed to set argument to the kernel! (ny)\n");
        return MHD_FAILURE;
    }
//	nz, z length
    status_in_setupCL = clSetKernelArg(
                fluidx_Kernel,
                4,
                sizeof(cl_int),
                (void *)&ctx_device.nz);
    if (status_in_setupCL != CL_SUCCESS)
    {
        printf("Error: Failed to set argument to the kernel! (nz)\n");
        return MHD_FAILURE;
    }
//	dt, cfl time 
    status_in_setupCL = clSetKernelArg(
                fluidx_Kernel,
                5,
                sizeof(cl_float),
                (void *)&ctx_device.dt);
    if (status_in_setupCL != CL_SUCCESS)
    {
        printf("Error: Failed to set argument to the kernel! (dt)\n");
        return MHD_FAILURE;
    }
//	adv_tmpB, matrix for b update
    status_in_setupCL = clSetKernelArg(
                fluidx_Kernel,
                6,
                sizeof(cl_mem),
                (void *)&ctx_device.adv_tmpB);
    if (status_in_setupCL != CL_SUCCESS)
    {
        printf("Error: Failed to set argument to the kernel! (adv_tmpB)\n");
        return MHD_FAILURE;
    }
//	u_update, fluid update
    status_in_setupCL = clSetKernelArg(
                fluidx_Kernel,
                7,
                sizeof(cl_mem),
                (void *)&ctx_device.u_update);
    if (status_in_setupCL != CL_SUCCESS)
    {
        printf("Error: Failed to set argument to the kernel! (u_update)\n");
        return MHD_FAILURE;
    }
//	b_update, magnetic update
    status_in_setupCL = clSetKernelArg(
                fluidx_Kernel,
                8,
                sizeof(cl_mem),
                (void *)&ctx_device.b_update);
    if (status_in_setupCL != CL_SUCCESS)
    {
        printf("Error: Failed to set argument to the kernel! (b_update)\n");
        return MHD_FAILURE;
    }
//	do the memory check
    status_in_setupCL = clGetKernelWorkGroupInfo(
		fluidx_Kernel,
		mhdComputeDevice[0],
                CL_KERNEL_LOCAL_MEM_SIZE,
                sizeof(cl_ulong),
                &mhdUsedLocalMemory,
                NULL);
printf("mhdUsedLocalMemory is %i\n",mhdUsedLocalMemory);
printf("mhdTotalLocalMemory is %i\n",mhdTotalLocalMemory);
    if (status_in_setupCL != CL_SUCCESS)
    {
	printf("Error: Failed to get CL_KERNEL_LOCAL_MEM_SIZE!\n");
        return MHD_FAILURE;
    }
    if(mhdUsedLocalMemory > mhdTotalLocalMemory)
    {
        printf("Unsupported: Insufficient local memory on device!\n");
        return MHD_FAILURE;
    }
//	get the other parameters about the kernel
    status_in_setupCL = clGetKernelWorkGroupInfo(
                fluidx_Kernel,
                mhdComputeDevice[0],
                CL_KERNEL_WORK_GROUP_SIZE,
                sizeof(size_t),
                &mhdKernelWorkGroupSize,
                NULL);
printf("mhdKernelWorkGroupSize is %i\n",mhdKernelWorkGroupSize);
    if (status_in_setupCL != CL_SUCCESS)
    {
        printf("Error: Failed to get CL_KERNEL_WORK_GROUP_SIZE!\n");
        return MHD_FAILURE;
    }
//
		return MHD_SUCCESS;
}

int
MHD::runCL_Kernel()
{
    cl_int   status_in_runCL = CL_SUCCESS;

    /*
 *      * Enqueue a kernel run call.
 *           */
    size_t globalThreads[] = {512};//{4*mhdKernelWorkGroupSize}; 
    size_t localThreads[] = {16};//{mhdKernelWorkGroupSize};

    if(localThreads[0] > mhdMaxWorkItemSizes[0] || localThreads[0] > mhdMaxWorkGroupSize)
    {
	printf("Unsupported: Device does not support requested number of work items!\n");
        return MHD_FAILURE;
    }

    status_in_runCL = clEnqueueNDRangeKernel(
		mhdCommandQueue,
                fluidx_Kernel,
                1,
                NULL,
                globalThreads,
                localThreads,
                0,
                NULL,
                NULL);
    if (status_in_runCL != CL_SUCCESS)
    {
	printf("Error: Failed to enqueue a kernel!\n");
	return MHD_FAILURE;
    }
//
    status_in_runCL = clFinish(mhdCommandQueue);
    if (status_in_runCL != CL_SUCCESS)
    {
	printf("Error: Failed to run clFinish!\n");
        return MHD_FAILURE;
    }

	return MHD_SUCCESS;
}

int 
MHD::setup()
{
//
    if(setupMHD_initialization() != MHD_SUCCESS)
	return MHD_FAILURE;
printf("finished setupMHD_initialization in setup()!\n");
//
    if(setupCL_ContextCommandMemory() != MHD_SUCCESS)
	return MHD_FAILURE;
printf("finished setupCL_ContextCommandMemory in setup()!\n");
//
    if(setupCL_ProgramKernel() != MHD_SUCCESS)
	return MHD_FAILURE;
printf("finished setupCL_ProgramKernel in setup()!\n");
//
    if(setupCL_ArgumentKernel() != MHD_SUCCESS)
	return MHD_FAILURE;
printf("finished setupCL_ArgumentKernel in setup()!\n");
//
	return MHD_SUCCESS;
}

int MHD::run()
{
//
    if(runCL_Kernel() != MHD_SUCCESS)
	return MHD_FAILURE;
printf("finished runCL_Kernel in run()!\n");
//
	return MHD_SUCCESS;
}

int
MHD::cleanup()
{
    /* Releases OpenCL resources (Context, Memory etc.) */
    cl_int status_in_cleanup = CL_SUCCESS;

    status_in_cleanup = clReleaseKernel(fluidx_Kernel);
    if (status_in_cleanup != CL_SUCCESS)
    {
	printf("Error: Failed to release kernel!\n");
	return MHD_FAILURE;
    }

    status_in_cleanup = clReleaseProgram(fluidx_Program);
    if (status_in_cleanup != CL_SUCCESS)
    {
        printf("Error: Failed to release program!\n");
        return MHD_FAILURE;
    }

    status_in_cleanup = clReleaseMemObject(ctx_device.u);
    if (status_in_cleanup != CL_SUCCESS)
    {
        printf("Error: Failed to release device memory! (ctx_device.u)\n");
        return MHD_FAILURE;
    }

    status_in_cleanup = clReleaseMemObject(ctx_device.b);
    if (status_in_cleanup != CL_SUCCESS)
    {
        printf("Error: Failed to release device memory! (ctx_device.b)\n");
        return MHD_FAILURE;
    }

    status_in_cleanup = clReleaseMemObject(ctx_device.adv_tmpB);
    if (status_in_cleanup != CL_SUCCESS)
    {
        printf("Error: Failed to release device memory! (ctx_device.adv_tmpB)\n");
        return MHD_FAILURE;
    }

    status_in_cleanup = clReleaseMemObject(ctx_device.u_update);
    if (status_in_cleanup != CL_SUCCESS)
    {
        printf("Error: Failed to release device memory! (ctx_device.u_update)\n");
        return MHD_FAILURE;
    }

    status_in_cleanup = clReleaseMemObject(ctx_device.b_update);
    if (status_in_cleanup != CL_SUCCESS)
    {
        printf("Error: Failed to release device memory! (ctx_device.b_update)\n");
        return MHD_FAILURE;
    }

    status_in_cleanup = clReleaseCommandQueue(mhdCommandQueue);
    if (status_in_cleanup != CL_SUCCESS)
    {
        printf("Error: Failed to release command queue!\n");
        return MHD_FAILURE;
    }

    status_in_cleanup = clReleaseContext(mhdComputeContext);
    if (status_in_cleanup != CL_SUCCESS)
    {
        printf("Error: Failed to release compute context!\n");
        return MHD_FAILURE;
    }
//
	return MHD_SUCCESS;
}


int main(int argc, char** argv)
{
printf("program start!\n");
// 
    MHD cl_MHD;
//
    if(cl_MHD.setup() != MHD_SUCCESS)
	return MHD_FAILURE;
//
    if(cl_MHD.run() != MHD_SUCCESS)
	return MHD_FAILURE;
//
    if (cl_MHD.cleanup() != MHD_SUCCESS)
	return MHD_FAILURE;
//
printf("program finish!\n");
	return EXIT_SUCCESS;
}
