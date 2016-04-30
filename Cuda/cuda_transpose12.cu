#include <stdio.h>
#include <math.h>
#include "cuda.h"

#include "parameter.h" 
#include "array_definition.h"
#include "cuda_funclist.h"

__global__ void transpose12_kernel(float *t12_ut, float *t12_bt, float *t12_u, float *t12_b, int *t12_nx, int *t12_ny, int *t12_nz);
//	

void cuda_transpose12(float *t12_ut, float *t12_bt, float *t12_u, float *t12_b, int *t12_nx, int *t12_ny, int *t12_nz, int *h_t12_nx, int *h_t12_ny, int *h_t12_nz)
{
//      initialization
int Totalthreads = (*h_t12_nx)*(*h_t12_ny)*(*h_t12_nz);
int numThreadsPerBlock = *h_t12_nx;
int numBlocks = Totalthreads/numThreadsPerBlock;
int NumOfU = 5;
int NumOfB = 3;
size_t u_memSize = NumOfU * numBlocks * numThreadsPerBlock * sizeof(float);
size_t b_memSize = NumOfB * numBlocks * numThreadsPerBlock * sizeof(float);
//      send it to device to calculate
dim3 dimGrid(*h_t12_ny,*h_t12_nz);
dim3 dimBlock(*h_t12_nx);
transpose12_kernel<<< dimGrid, dimBlock >>>(t12_ut,t12_bt,t12_u,t12_b,t12_nx,t12_ny,t12_nz);
//
cudaThreadSynchronize();
//
checkCUDAError("kernel execution in cuda_transpose12");
//	cudaMemcpy
//	from d_ut to d_u, in device
//	from d_bt to d_b, in device
cudaMemcpy(t12_u,t12_ut, u_memSize, cudaMemcpyDeviceToDevice );
cudaMemcpy(t12_b,t12_bt, b_memSize, cudaMemcpyDeviceToDevice );
//
checkCUDAError("memcpy: from device to device, in cuda_transpose12");
//
}

__global__ void transpose12_kernel(float *t12_ut, float *t12_bt, float *t12_u, float *t12_b, int *t12_nx, int *t12_ny, int *t12_nz)
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
t12_ut[a4D_FinC(5,gridDim.x,blockDim.x,gridDim.y,(1-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t12_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,(1-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t12_ut[a4D_FinC(5,gridDim.x,blockDim.x,gridDim.y,(2-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t12_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,(3-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t12_ut[a4D_FinC(5,gridDim.x,blockDim.x,gridDim.y,(3-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t12_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,(2-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t12_ut[a4D_FinC(5,gridDim.x,blockDim.x,gridDim.y,(4-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t12_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,(4-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t12_ut[a4D_FinC(5,gridDim.x,blockDim.x,gridDim.y,(5-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t12_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,(5-1),threadIdx.x,blockIdx.x,blockIdx.y)];
//
t12_bt[a4D_FinC(3,gridDim.x,blockDim.x,gridDim.y,(1-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t12_b[a4D_FinC(3,blockDim.x,gridDim.x,gridDim.y,(2-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t12_bt[a4D_FinC(3,gridDim.x,blockDim.x,gridDim.y,(2-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t12_b[a4D_FinC(3,blockDim.x,gridDim.x,gridDim.y,(1-1),threadIdx.x,blockIdx.x,blockIdx.y)];
t12_bt[a4D_FinC(3,gridDim.x,blockDim.x,gridDim.y,(3-1),blockIdx.x,threadIdx.x,blockIdx.y)]=t12_b[a4D_FinC(3,blockDim.x,gridDim.x,gridDim.y,(3-1),threadIdx.x,blockIdx.x,blockIdx.y)];
return;
}

//	
//-----------------------
/*
   Fortran subroutine arguments are passed by references.
   call fun( array_a, array_b, N) will be mapped to
   function (*a, *b, *N);
*/
extern "C" void cuda_transpose12_(float *h_ut, float *h_bt, float *h_u, float *h_b, int *h_nx, int *h_ny, int *h_nz, float *h_dt)
{
int Totalthreads = (*h_nx)*(*h_ny)*(*h_nz);
int numThreadsPerBlock = *h_nx;
int numBlocks = Totalthreads/numThreadsPerBlock;
int NumOfU = 5;
int NumOfB = 3;
//      intialize
size_t u_memSize = NumOfU * numBlocks * numThreadsPerBlock * sizeof(float);
size_t b_memSize = NumOfB * numBlocks * numThreadsPerBlock * sizeof(float);
//
float *d_u, *d_b;
cudaMalloc( (void **) &d_u, u_memSize );
cudaMalloc( (void **) &d_b, b_memSize );
cudaMemcpy( d_u, h_u, u_memSize, cudaMemcpyHostToDevice );
cudaMemcpy( d_b, h_b, b_memSize, cudaMemcpyHostToDevice );
//
float *d_ut, *d_bt;
cudaMalloc( (void **) &d_ut, u_memSize );
cudaMalloc( (void **) &d_bt, b_memSize );
//
int *d_nx,*d_ny,*d_nz;
size_t n_memSize = sizeof(int);
cudaMalloc( (void **) &d_nx, n_memSize );
cudaMalloc( (void **) &d_ny, n_memSize );
cudaMalloc( (void **) &d_nz, n_memSize );
cudaMemcpy( d_nx, h_nx, n_memSize, cudaMemcpyHostToDevice );
cudaMemcpy( d_ny, h_ny, n_memSize, cudaMemcpyHostToDevice );
cudaMemcpy( d_nz, h_nz, n_memSize, cudaMemcpyHostToDevice );
//
dim3 dimGrid(*h_ny,*h_nz);
dim3 dimBlock(numThreadsPerBlock);
transpose12_kernel<<< dimGrid, dimBlock >>>( d_ut, d_bt, d_u, d_b, d_nx, d_ny, d_nz);
//
cudaThreadSynchronize();
//
checkCUDAError("kernel execution");
//
//	find the max
cudaMemcpy( h_ut, d_ut, u_memSize, cudaMemcpyDeviceToHost );
cudaMemcpy( h_bt, d_bt, b_memSize, cudaMemcpyDeviceToHost );
//
checkCUDAError("memcpy");
//
cudaFree(d_u);
cudaFree(d_b);
cudaFree(d_nx);
cudaFree(d_ny);
cudaFree(d_nz);
cudaFree(d_ut);
cudaFree(d_bt);
//
return;
}



