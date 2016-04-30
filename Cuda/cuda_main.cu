#include <stdio.h>
#include <math.h>
#include "cuda.h"

#include "parameter.h" 
#include "array_definition.h"
#include "cuda_funclist.h"
#include "cuda_function.h" 
#include "cuda_subroutine.h"

extern "C" void cuda_main(float *h_u, float *h_b, int *h_nx, int *h_ny, int *h_nz)
{
//	general info initialization
int Totalthreads = (*h_nx)*(*h_ny)*(*h_nz);
int numThreadsPerBlock = *h_nx;
int numBlocks = Totalthreads/numThreadsPerBlock;
int NumOfU = 5;
int NumOfB = 3;
//	memory size initialization
size_t u_memSize = NumOfU * numBlocks * numThreadsPerBlock * sizeof(float);
size_t b_memSize = NumOfB * numBlocks * numThreadsPerBlock * sizeof(float);
size_t c_memSize = numBlocks * numThreadsPerBlock * sizeof(float);
size_t int_memSize = sizeof(int);
size_t float_memSize = sizeof(float);
//	data on the host
float *h_dt;
//	data on the device
//	cudaMalloc
//	for general purpose
float *d_u, *d_b;
cudaMalloc( (void **) &d_u, u_memSize );
cudaMalloc( (void **) &d_b, b_memSize );
int *d_nx,*d_ny,*d_nz;
cudaMalloc( (void **) &d_nx, int_memSize );
cudaMalloc( (void **) &d_ny, int_memSize );
cudaMalloc( (void **) &d_nz, int_memSize );
float *d_dt;
cudaMalloc( (void **) &d_dt, float_memSize );
//	for cuda_cfl
float *d_c;
cudaMalloc( (void **) &d_c, c_memSize );
//	for cuda_advectbyzx
float *d_temp;
cudaMalloc( (void **) &d_temp, c_memSize );
//	for cuda_transpose
float *d_ut, *d_bt;
cudaMalloc( (void **) &d_ut, u_memSize );
cudaMalloc( (void **) &d_bt, b_memSize );
//	cudaMemcpy
//	copy data from host to device
cudaMemcpy( d_u, h_u, u_memSize, cudaMemcpyHostToDevice );
cudaMemcpy( d_b, h_b, b_memSize, cudaMemcpyHostToDevice );
cudaMemcpy( d_nx, h_nx, int_memSize, cudaMemcpyHostToDevice );
cudaMemcpy( d_ny, h_ny, int_memSize, cudaMemcpyHostToDevice );
cudaMemcpy( d_nz, h_nz, int_memSize, cudaMemcpyHostToDevice );
//
checkCUDAError("memcpy: from host to device, in cuda_main");
//	initialize data for loop
float t,dt,tf;
int iter;
float ct;
t=0;
iter=0;
ct=100.;
tf=ct*10;
//	initialization for timing
cudaEvent_t start, stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);
//	in milliseconds with a resolution of around 0.5 microseconds
float elapsedTime;
do {
//	start the timer
	cudaEventRecord(start,0);
//	output
//	if you want to output, you have to use cudaMemcpy
//	copy the data from device to host to output
	cudaMemcpy( h_u, d_u, u_memSize, cudaMemcpyDeviceToHost );
	cudaMemcpy( h_b, d_b, b_memSize, cudaMemcpyDeviceToHost );
	printf("t=	%f,	%i,	%f\n",t,iter,h_u[a4D_FinC(5,(*h_nx),(*h_ny),(*h_nz),(5-1),(*h_nx)/4,1,1)]);
//	done output
	iter=iter+1;
	cuda_cfl(d_u,d_b,d_nx,d_ny,d_nz,d_dt,d_c,h_nx,h_ny,h_nz,h_dt);
	dt=0.9*(*h_dt);
	//dt=0.5;
	if (dt>(tf-t)/2.0) dt=(tf-t)/2.0;
	t=t+2.0*dt;
//	start sweep
	cuda_fluidx(d_u,d_b,d_nx,d_ny,d_nz,d_dt,h_nx,h_ny,h_nz);
	cuda_advectbyzx(d_u,d_b,d_nx,d_ny,d_nz,d_dt,d_temp,h_nx,h_ny,h_nz);
//	the y sweep
	cuda_transpose12(d_ut,d_bt,d_u,d_b,d_nx,d_ny,d_nz,h_nx,h_ny,h_nz);
	cuda_fluidx(d_u,d_b,d_ny,d_nx,d_nz,d_dt,h_ny,h_nx,h_nz);
	cuda_advectbyzx(d_u,d_b,d_ny,d_nx,d_nz,d_dt,d_temp,h_ny,h_nx,h_nz);
//	z sweep
	cuda_transpose13(d_ut,d_bt,d_u,d_b,d_ny,d_nx,d_nz,h_ny,h_nx,h_nz);
	cuda_fluidx(d_u,d_b,d_nz,d_nx,d_ny,d_dt,h_nz,h_nx,h_ny);
	cuda_advectbyzx(d_u,d_b,d_nz,d_nx,d_ny,d_dt,d_temp,h_nz,h_nx,h_ny);
	cuda_advectbyzx(d_u,d_b,d_nz,d_nx,d_ny,d_dt,d_temp,h_nz,h_nx,h_ny);
	cuda_fluidx(d_u,d_b,d_nz,d_nx,d_ny,d_dt,h_nz,h_nx,h_ny);

//	back
	cuda_transpose13(d_ut,d_bt,d_u,d_b,d_nz,d_nx,d_ny,h_nz,h_nx,h_ny);
	cuda_advectbyzx(d_u,d_b,d_ny,d_nx,d_nz,d_dt,d_temp,h_ny,h_nx,h_nz);
	cuda_fluidx(d_u,d_b,d_ny,d_nx,d_nz,d_dt,h_ny,h_nx,h_nz);
//	x again
	cuda_transpose12(d_ut,d_bt,d_u,d_b,d_ny,d_nx,d_nz,h_ny,h_nx,h_nz);
	cuda_advectbyzx(d_u,d_b,d_nx,d_ny,d_nz,d_dt,d_temp,h_nx,h_ny,h_nz);
	cuda_fluidx(d_u,d_b,d_nx,d_ny,d_nz,d_dt,h_nx,h_ny,h_nz);
//	finish sweep
//	stop the timer
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime,start,stop);	
	printf("time per loop(in milliseconds):	%f\n",elapsedTime);
} while (t<tf);
//
//      cudaMemcpy
//      copy data from device to host
cudaMemcpy( h_u, d_u, u_memSize, cudaMemcpyDeviceToHost );
cudaMemcpy( h_b, d_b, b_memSize, cudaMemcpyDeviceToHost );
//
checkCUDAError("memcpy: from device to host, in cuda_main");
//
cudaFree(d_u);
cudaFree(d_b);
cudaFree(d_nx);
cudaFree(d_ny);
cudaFree(d_nz);
cudaFree(d_dt);
//
cudaEventDestroy(start);
cudaEventDestroy(stop);
//
return;
}

