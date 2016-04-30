#include <stdio.h>
#include <math.h>
#include "cuda.h"

#include "parameter.h" 
#include "array_definition.h"
#include "cuda_funclist.h"

__global__ void cfl_kernel(float *cfl_u, float *cfl_b, int *cfl_nx, int *cfl_ny, int *cfl_nz, float *cfl_dt, float *cfl_c);
__device__ float max_3num(float *m3_1, float *m3_2, float *m3_3);
__device__ float max_2num(float *m2_1, float *m2_2);
__host__ void h_cfl_find_max(float *hcfm_in, int *hcfm_ny, int *hcfm_nz, float *hcfm_out);

void cuda_cfl(float *cfl_u, float *cfl_b, int *cfl_nx, int *cfl_ny, int *cfl_nz, float *cfl_dt, float *cfl_c, int *h_cfl_nx, int *h_cfl_ny, int *h_cfl_nz, float *h_cfl_dt)
{
//	initialization
int Totalthreads = (*h_cfl_nx)*(*h_cfl_ny)*(*h_cfl_nz);
int numThreadsPerBlock = *h_cfl_nx;
int numBlocks = Totalthreads/numThreadsPerBlock;
size_t c_memSize = numBlocks * sizeof(float);
size_t dt_memSize = sizeof(float);
//	send it to device to calculate
dim3 dimGrid(*h_cfl_ny,*h_cfl_nz);
dim3 dimBlock(*h_cfl_nx);
cfl_kernel<<< dimGrid, dimBlock >>>( cfl_u, cfl_b, cfl_nx, cfl_ny, cfl_nz, cfl_dt, cfl_c);
//
cudaThreadSynchronize();
//
checkCUDAError("kernel execution in cuda_cfl");
//      get it from device to find the max of c
float *temp_h_cfl_c;
temp_h_cfl_c = (float *) malloc(c_memSize);
cudaMemcpy( temp_h_cfl_c, cfl_c, c_memSize, cudaMemcpyDeviceToHost );
//
checkCUDAError("memcpy: from device to host, in cuda_cfl");
//
float max_c;
max_c=0;
h_cfl_find_max(temp_h_cfl_c,h_cfl_ny,h_cfl_nz,&max_c);
//	find it and get cfl_dt in host
(*h_cfl_dt)=1/(max_c);
//	copy it to device
cudaMemcpy( cfl_dt, h_cfl_dt, dt_memSize, cudaMemcpyHostToDevice );
//
checkCUDAError("memcpy: from host to device, in cuda_cfl");
//
free(temp_h_cfl_c);
//
}

__global__ void cfl_kernel(float *cfl_u, float *cfl_b, int *cfl_nx, int *cfl_ny, int *cfl_nz, float *cfl_dt, float *cfl_c)
{
/*
two dimensional array of blocks on grid where each block has one dimensional array of threads:
UniqueBlockIndex = blockIdx.y * gridDim.x + blockIdx.x;
UniqueThreadIndex = UniqueBlockIndex * blockDim.x + threadIdx.x;
*/
__shared__ float cfl_s_b1_ip[BLOCK_SIZE];
__shared__ float cfl_s_b2_jp[BLOCK_SIZE];
__shared__ float cfl_s_b3_kp[BLOCK_SIZE];
__shared__ float cfl_s_b[3*BLOCK_SIZE];
__shared__ float cfl_s_u[5*BLOCK_SIZE];
/*
i = threadIdx.x
j = blockIdx.x
k = blockIdx.y
nx = blockDim.x
ny = gridDim.x
nz = gridDim.y
*/
float gamma;
gamma=5.0/3.0;
int ii;
int kp,jp,ip;
//      kp=mod(k,nz)+1
kp=(blockIdx.y+1)%(*cfl_nz);
//      jp=mod(j,ny)+1
jp=(blockIdx.x+1)%(*cfl_ny);
//      ip=mod(i,nx)+1
ip=(threadIdx.x+1)%(*cfl_nx);
//	get cfl_s_b1_ip
cfl_s_b1_ip[threadIdx.x]=cfl_b[a4D_FinC(3,blockDim.x,gridDim.x,gridDim.y,(1-1),ip,blockIdx.x,blockIdx.y)];
//	get cfl_s_b2_jp
cfl_s_b2_jp[threadIdx.x]=cfl_b[a4D_FinC(3,blockDim.x,gridDim.x,gridDim.y,(2-1),threadIdx.x,jp,blockIdx.y)];
//	get cfl_s_b3_kp
cfl_s_b3_kp[threadIdx.x]=cfl_b[a4D_FinC(3,blockDim.x,gridDim.x,gridDim.y,(3-1),threadIdx.x,blockIdx.x,kp)];
//	get cfl_s_u
for (ii=0;ii<5;ii++)
{
	cfl_s_u[a2D_FinC(5,blockDim.x,ii,threadIdx.x)]=cfl_u[a4D_FinC(5,blockDim.x,gridDim.x,gridDim.y,ii,threadIdx.x,blockIdx.x,blockIdx.y)];
}
//	get cfl_s_b
for (ii=0;ii<3;ii++)
{
	cfl_s_b[a2D_FinC(3,blockDim.x,ii,threadIdx.x)]=cfl_b[a4D_FinC(3,blockDim.x,gridDim.x,gridDim.y,ii,threadIdx.x,blockIdx.x,blockIdx.y)];
}
//
__syncthreads();
//
float bx,by,bz;
//	bx=(b(1,i,j,k)+b(1,ip,j,k))/2
bx=(cfl_s_b[a2D_FinC(3,blockDim.x,(1-1),threadIdx.x)]+cfl_s_b1_ip[threadIdx.x])/2.0;
//	by=(b(2,i,j,k)+b(2,i,jp,k))/2
by=(cfl_s_b[a2D_FinC(3,blockDim.x,(2-1),threadIdx.x)]+cfl_s_b2_jp[threadIdx.x])/2.0;
//	bz=(b(3,i,j,k)+b(3,i,j,kp))/2
bz=(cfl_s_b[a2D_FinC(3,blockDim.x,(3-1),threadIdx.x)]+cfl_s_b3_kp[threadIdx.x])/2.0;
float v;
//	v=maxval(abs(u(2:4,i,j,k)/u(1,i,j,k)))
float temp1,temp2,temp3;
temp1=fabs(cfl_s_u[a2D_FinC(5,blockDim.x,(2-1),threadIdx.x)]/cfl_s_u[a2D_FinC(5,blockDim.x,(1-1),threadIdx.x)]);
temp2=fabs(cfl_s_u[a2D_FinC(5,blockDim.x,(3-1),threadIdx.x)]/cfl_s_u[a2D_FinC(5,blockDim.x,(1-1),threadIdx.x)]);
temp3=fabs(cfl_s_u[a2D_FinC(5,blockDim.x,(4-1),threadIdx.x)]/cfl_s_u[a2D_FinC(5,blockDim.x,(1-1),threadIdx.x)]);
v=max_3num(&temp1,&temp2,&temp3);
float b2;
b2=bx*bx+by*by+bz*bz;
//	ps=(u(5,i,j,k)-sum(u(2:4,i,j,k)**2,1)/u(1,i,j,k)/2)*(gamma-1)+(2-gamma)*b2/2
float ps;
ps=(cfl_s_u[a2D_FinC(5,blockDim.x,(5-1),threadIdx.x)]-(cfl_s_u[a2D_FinC(5,blockDim.x,(2-1),threadIdx.x)]*cfl_s_u[a2D_FinC(5,blockDim.x,(2-1),threadIdx.x)]+cfl_s_u[a2D_FinC(5,blockDim.x,(3-1),threadIdx.x)]*cfl_s_u[a2D_FinC(5,blockDim.x,(3-1),threadIdx.x)]+cfl_s_u[a2D_FinC(5,blockDim.x,(4-1),threadIdx.x)]*cfl_s_u[a2D_FinC(5,blockDim.x,(4-1),threadIdx.x)])/cfl_s_u[a2D_FinC(5,blockDim.x,(1-1),threadIdx.x)]/2.0)*(gamma-1.0)+(2.0-gamma)*b2/2.0;
//	p=ps-b2/2
float p;
p=ps-b2/2.0;
//	c=max(c,v+sqrt(abs(  (b2*2+gamma*p)/u(1,i,j,k))))
temp1=v+sqrt(fabs((b2*2.0+gamma*p)/cfl_s_u[a2D_FinC(5,blockDim.x,(1-1),threadIdx.x)]));
//	find max
__shared__ float cfl_s_c[BLOCK_SIZE];
float temp_c_max;
cfl_s_c[threadIdx.x]=temp1;
__syncthreads();
if (threadIdx.x==0)
{
temp_c_max=0.0;
for (int i=0; i<BLOCK_SIZE; i++)
{
        if (cfl_s_c[i]>temp_c_max) temp_c_max=cfl_s_c[i];
}
(*cfl_c)=temp_c_max;
}
//
return;
}

__device__ float max_3num(float *m3_1, float *m3_2, float *m3_3)
{
if ((*m3_1)>(*m3_2))
{
        if ((*m3_1)>(*m3_3))
        {
                return (*m3_1);
        }
        else
        {
                return (*m3_3);
        }
}
else
{
        if ((*m3_2)>(*m3_3))
        {
                return (*m3_2);
        }
        else
        {
                return (*m3_3);
        }
}
}

__device__ float max_2num(float *m2_1, float *m2_2)
{
if ((*m2_1)>(*m2_2))
{
        return (*m2_1);
}
else
{
        return (*m2_2);
}
}

__host__ void h_cfl_find_max(float *hcfm_in, int *hcfm_ny, int *hcfm_nz, float *hcfm_out)
{
int j,k;
(*hcfm_out)=0;
for (k=0;k<(*hcfm_nz);k++)
{
	for (j=0;j<(*hcfm_ny);j++)
	{
		if (hcfm_in[a2D_FinC((*hcfm_ny),(*hcfm_nz),j,k)]>(*hcfm_out))
		{
			(*hcfm_out)=hcfm_in[a2D_FinC((*hcfm_ny),(*hcfm_nz),j,k)];
		}		
	}
}
}

