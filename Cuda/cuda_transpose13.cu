#include <stdio.h>
#include <math.h>
#include "cuda.h"

#include "parameter.h" 
#include "array_definition.h"
#include "cuda_funclist.h"

__global__ void transpose13_kernel(float *t13_ut, float *t13_bt, float *t13_u, float *t13_b, int *t13_nx, int *t13_ny, int *t13_nz);
__global__ void transpose12_kernel(float *t12_ut, float *t12_bt, float *t12_u, float *t12_b, int *t12_nx, int *t12_ny, int *t12_nz);
__global__ void transpose13_nz1_kernel(float *t13_nz1_ut, float *t13_nz1_bt, float *t13_nz1_u, float *t13_nz1_b, int *t13_nz1_nx, int *t13_nz1_ny, int *t13_nz1_nz);
__global__ void transpose13_nx1_kernel(float *t13_nx1_ut, float *t13_nx1_bt, float *t13_nx1_u, float *t13_nx1_b, int *t13_nx1_nx, int *t13_nx1_ny, int *t13_nx1_nz);

void cuda_transpose13(float *t13_ut, float *t13_bt, float *t13_u, float *t13_b, int *t13_nx, int *t13_ny, int *t13_nz, int *h_t13_nx, int *h_t13_ny, int *h_t13_nz)
{
//      initialization
int Totalthreads = (*h_t13_nx)*(*h_t13_ny)*(*h_t13_nz);
int numThreadsPerBlock = *h_t13_nx;
int numBlocks = Totalthreads/numThreadsPerBlock;
int NumOfU = 5;
int NumOfB = 3;
size_t u_memSize = NumOfU * numBlocks * numThreadsPerBlock * sizeof(float);
size_t b_memSize = NumOfB * numBlocks * numThreadsPerBlock * sizeof(float);
//      send it to device to calculate
dim3 dimGrid(*h_t13_ny,*h_t13_nz);
dim3 dimBlock(*h_t13_nx);
if ((*h_t13_nx)==(*h_t13_nz))
{
        transpose13_kernel<<< dimGrid, dimBlock >>>(t13_ut,t13_bt,t13_u,t13_b,t13_nx,t13_ny,t13_nz);
}
else if ((*h_t13_nz)==1)
{
        transpose13_nz1_kernel<<< dimGrid, dimBlock >>>(t13_ut,t13_bt,t13_u,t13_b,t13_nx,t13_ny,t13_nz);
}
else if ((*h_t13_nx)==1)
{
        transpose13_nx1_kernel<<< dimGrid, dimBlock >>>(t13_ut,t13_bt,t13_u,t13_b,t13_nx,t13_ny,t13_nz);
}
else
{
        printf("nz<>nx not supported\n");
}
//
cudaThreadSynchronize();
//
checkCUDAError("kernel execution in cuda_transpose13");
//      cudaMemcpy
//      from d_ut to d_u, in device
//      from d_bt to d_b, in device
cudaMemcpy(t13_u,t13_ut, u_memSize, cudaMemcpyDeviceToDevice );
cudaMemcpy(t13_b,t13_bt, b_memSize, cudaMemcpyDeviceToDevice );
//
checkCUDAError("memcpy: from device to device, in cuda_transpose13");
//
}

__global__ void transpose13_kernel(float *t13_ut, float *t13_bt, float *t13_u, float *t13_b, int *t13_nx, int *t13_ny, int *t13_nz)
{
/*
two dimensional array of blocks on grid where each block has one dimensional array of threads:
UniqueBlockIndex = blockIdx.y * gridDim.x + blockIdx.x;
UniqueThreadIndex = UniqueBlockIndex * blockDim.x + threadIdx.x;
*/
/*
i = threadIdx.x
j = blockIdx.x
k = blockIdx.y
nx = blockDim.x
ny = gridDim.x
nz = gridDim.y
*/
t13_ut[a4D_FinC(5,gridDim.y,gridDim.x,blockDim.x,(1-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,(1-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t13_ut[a4D_FinC(5,gridDim.y,gridDim.x,blockDim.x,(2-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,(4-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t13_ut[a4D_FinC(5,gridDim.y,gridDim.x,blockDim.x,(3-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,(3-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t13_ut[a4D_FinC(5,gridDim.y,gridDim.x,blockDim.x,(4-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,(2-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t13_ut[a4D_FinC(5,gridDim.y,gridDim.x,blockDim.x,(5-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,(5-1),threadIdx.x,blockIdx.x,blockIdx.y)];
//
t13_bt[a4D_FinC(3,gridDim.y,gridDim.x,blockDim.x,(1-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_b[a4D_FinC(3,blockDim.x,gridDim.x,gridDim.y,(3-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t13_bt[a4D_FinC(3,gridDim.y,gridDim.x,blockDim.x,(2-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_b[a4D_FinC(3,blockDim.x,gridDim.x,gridDim.y,(2-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t13_bt[a4D_FinC(3,gridDim.y,gridDim.x,blockDim.x,(3-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_b[a4D_FinC(3,blockDim.x,gridDim.x,gridDim.y,(1-1),threadIdx.x,blockIdx.x,blockIdx.y)];
return;
}

__global__ void transpose13_nz1_kernel(float *t13_nz1_ut, float *t13_nz1_bt, float *t13_nz1_u, float *t13_nz1_b, int *t13_nz1_nx, int *t13_nz1_ny, int *t13_nz1_nz)
{
/*
i = threadIdx.x
j = blockIdx.x
k = blockIdx.y
nx = blockDim.x
ny = gridDim.x
nz = gridDim.y
*/
//	transpose12
t13_nz1_ut[a4D_FinC(5,gridDim.x,blockDim.x,gridDim.y,(1-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t13_nz1_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,(1-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t13_nz1_ut[a4D_FinC(5,gridDim.x,blockDim.x,gridDim.y,(2-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t13_nz1_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,(3-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t13_nz1_ut[a4D_FinC(5,gridDim.x,blockDim.x,gridDim.y,(3-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t13_nz1_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,(2-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t13_nz1_ut[a4D_FinC(5,gridDim.x,blockDim.x,gridDim.y,(4-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t13_nz1_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,(4-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t13_nz1_ut[a4D_FinC(5,gridDim.x,blockDim.x,gridDim.y,(5-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t13_nz1_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,(5-1),threadIdx.x,blockIdx.x,blockIdx.y)];
//
t13_nz1_bt[a4D_FinC(3,gridDim.x,blockDim.x,gridDim.y,(1-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t13_nz1_b[a4D_FinC(3,blockDim.x,gridDim.x,gridDim.y,(2-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t13_nz1_bt[a4D_FinC(3,gridDim.x,blockDim.x,gridDim.y,(2-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t13_nz1_b[a4D_FinC(3,blockDim.x,gridDim.x,gridDim.y,(1-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t13_nz1_bt[a4D_FinC(3,gridDim.x,blockDim.x,gridDim.y,(3-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t13_nz1_b[a4D_FinC(3,blockDim.x,gridDim.x,gridDim.y,(3-1),threadIdx.x,blockIdx.x,blockIdx.y)];
//
//	second part
float temp1,temp2,temp3;
temp1=t13_nz1_ut[a4D_FinC(5,1,gridDim.x,blockDim.x,(2-1),(1-1),blockIdx.x,threadIdx.x)];
temp2=t13_nz1_ut[a4D_FinC(5,1,gridDim.x,blockDim.x,(3-1),(1-1),blockIdx.x,threadIdx.x)];
temp3=t13_nz1_ut[a4D_FinC(5,1,gridDim.x,blockDim.x,(4-1),(1-1),blockIdx.x,threadIdx.x)];
t13_nz1_ut[a4D_FinC(5,1,gridDim.x,blockDim.x,(2-1),(1-1),blockIdx.x,threadIdx.x)]=temp3;
t13_nz1_ut[a4D_FinC(5,1,gridDim.x,blockDim.x,(3-1),(1-1),blockIdx.x,threadIdx.x)]=temp1;
t13_nz1_ut[a4D_FinC(5,1,gridDim.x,blockDim.x,(4-1),(1-1),blockIdx.x,threadIdx.x)]=temp2;
//
temp1=t13_nz1_bt[a4D_FinC(3,1,gridDim.x,blockDim.x,(1-1),(1-1),blockIdx.x,threadIdx.x)];
temp2=t13_nz1_bt[a4D_FinC(3,1,gridDim.x,blockDim.x,(2-1),(1-1),blockIdx.x,threadIdx.x)];
temp3=t13_nz1_bt[a4D_FinC(3,1,gridDim.x,blockDim.x,(3-1),(1-1),blockIdx.x,threadIdx.x)];
t13_nz1_bt[a4D_FinC(3,1,gridDim.x,blockDim.x,(1-1),(1-1),blockIdx.x,threadIdx.x)]=temp3;
t13_nz1_bt[a4D_FinC(3,1,gridDim.x,blockDim.x,(1-1),(2-1),blockIdx.x,threadIdx.x)]=temp1;
t13_nz1_bt[a4D_FinC(3,1,gridDim.x,blockDim.x,(1-1),(3-1),blockIdx.x,threadIdx.x)]=temp2;
//
return;
}

__global__ void transpose13_nx1_kernel(float *t13_nx1_ut, float *t13_nx1_bt, float *t13_nx1_u, float *t13_nx1_b, int *t13_nx1_nx, int *t13_nx1_ny, int *t13_nx1_nz)
{
/*
i = threadIdx.x
j = blockIdx.x
k = blockIdx.y
nx = blockDim.x
ny = gridDim.x
nz = gridDim.y
*/
//      transpose12
t13_nx1_ut[a4D_FinC(5,gridDim.y,gridDim.x,blockDim.x,(1-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_nx1_u[a4D_FinC(5,gridDim.x,gridDim.y,blockDim.x,(1-1),blockIdx.x,blockIdx.y,threadIdx.x)];
t13_nx1_ut[a4D_FinC(5,gridDim.y,gridDim.x,blockDim.x,(2-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_nx1_u[a4D_FinC(5,gridDim.x,gridDim.y,blockDim.x,(3-1),blockIdx.x,blockIdx.y,threadIdx.x)];
t13_nx1_ut[a4D_FinC(5,gridDim.y,gridDim.x,blockDim.x,(3-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_nx1_u[a4D_FinC(5,gridDim.x,gridDim.y,blockDim.x,(2-1),blockIdx.x,blockIdx.y,threadIdx.x)];
t13_nx1_ut[a4D_FinC(5,gridDim.y,gridDim.x,blockDim.x,(4-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_nx1_u[a4D_FinC(5,gridDim.x,gridDim.y,blockDim.x,(4-1),blockIdx.x,blockIdx.y,threadIdx.x)];
t13_nx1_ut[a4D_FinC(5,gridDim.y,gridDim.x,blockDim.x,(5-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_nx1_u[a4D_FinC(5,gridDim.x,gridDim.y,blockDim.x,(5-1),blockIdx.x,blockIdx.y,threadIdx.x)];
//
t13_nx1_bt[a4D_FinC(3,gridDim.y,gridDim.x,blockDim.x,(1-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_nx1_b[a4D_FinC(3,gridDim.x,gridDim.y,blockDim.x,(2-1),blockIdx.x,blockIdx.y,threadIdx.x)];
t13_nx1_bt[a4D_FinC(3,gridDim.y,gridDim.x,blockDim.x,(2-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_nx1_b[a4D_FinC(3,gridDim.x,gridDim.y,blockDim.x,(1-1),blockIdx.x,blockIdx.y,threadIdx.x)];
t13_nx1_bt[a4D_FinC(3,gridDim.y,gridDim.x,blockDim.x,(3-1),blockIdx.y,blockIdx.x,threadIdx.x)]=t13_nx1_b[a4D_FinC(3,gridDim.x,gridDim.y,blockDim.x,(3-1),blockIdx.x,blockIdx.y,threadIdx.x)];
//
//	second part
float temp1,temp2,temp3;
temp1=t13_nx1_ut[a4D_FinC(5,gridDim.y,gridDim.x,1,(2-1),blockIdx.y,blockIdx.x,(1-1))];
temp2=t13_nx1_ut[a4D_FinC(5,gridDim.y,gridDim.x,1,(3-1),blockIdx.y,blockIdx.x,(1-1))];
temp3=t13_nx1_ut[a4D_FinC(5,gridDim.y,gridDim.x,1,(4-1),blockIdx.y,blockIdx.x,(1-1))];
t13_nx1_ut[a4D_FinC(5,gridDim.y,gridDim.x,1,(2-1),blockIdx.y,blockIdx.x,(1-1))]=temp3;
t13_nx1_ut[a4D_FinC(5,gridDim.y,gridDim.x,1,(3-1),blockIdx.y,blockIdx.x,(1-1))]=temp1;
t13_nx1_ut[a4D_FinC(5,gridDim.y,gridDim.x,1,(4-1),blockIdx.y,blockIdx.x,(1-1))]=temp2;
//
temp1=t13_nx1_bt[a4D_FinC(3,gridDim.y,gridDim.x,1,(1-1),blockIdx.y,blockIdx.x,(1-1))];
temp2=t13_nx1_bt[a4D_FinC(3,gridDim.y,gridDim.x,1,(2-1),blockIdx.y,blockIdx.x,(1-1))];
temp3=t13_nx1_bt[a4D_FinC(3,gridDim.y,gridDim.x,1,(3-1),blockIdx.y,blockIdx.x,(1-1))];
t13_nx1_bt[a4D_FinC(3,gridDim.y,gridDim.x,1,(1-1),blockIdx.y,blockIdx.x,(1-1))]=temp3;
t13_nx1_bt[a4D_FinC(3,gridDim.y,gridDim.x,1,(2-1),blockIdx.y,blockIdx.x,(1-1))]=temp1;
t13_nx1_bt[a4D_FinC(3,gridDim.y,gridDim.x,1,(3-1),blockIdx.y,blockIdx.x,(1-1))]=temp2;
//
return;
}
