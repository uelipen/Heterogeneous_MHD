#include "data_type.hpp"

__kernel 
void 
fluidx( 
	__global data_type_t* ctx_device.u,
	__global data_type_t* ctx_device.b,
	int ctx_device.nx,
	int ctx_device.ny,
	int ctx_device.nz,
	data_type_t ctx_device.dt,
	__global data_type_t* ctx_device.adv_tmpB,
	__global data_type_t* ctx_device.u_update,
	__global data_type_t* ctx_device.b_update)
{                                                                      

};

