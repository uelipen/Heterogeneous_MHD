#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#include <stdio.h>
#include <stdlib.h>

#include "data_type.h"

#include "array_definition.h"
#include "value_spu.h"

#include "funclist_spu.h"

void tvdb(data_type_t *tvdb_flux, data_type_t *tvdb_b, data_type_t *tvdb_vg, int *tvdb_n, data_type_t *tvdb_dt)
{

int i; 

data_type_t vh[*tvdb_n] __attribute__ ((aligned (16)));
data_type_t flux1[*tvdb_n] __attribute__ ((aligned (16)));
data_type_t b1[*tvdb_n] __attribute__ ((aligned (16)));

data_type_t tvdb_tmp1[*tvdb_n] __attribute__ ((aligned (16))); 
data_type_t tvdb_tmp2[*tvdb_n] __attribute__ ((aligned (16)));

int ip,ipp,im; 
data_type_t v,w,wp,wm,dw; 

//	vh=(vg+cshift(vg,1))/2
cshift_matrix1D(tvdb_vg,(*tvdb_n),1,tvdb_tmp1);
abXYc_matrix1D(vh,tvdb_vg,tvdb_tmp1,0.5E0,0.5E0,(*tvdb_n),0.0E0);
//	  where(vh>0)
//	     flux1=b*vg
//	  elsewhere
for (i=0;i<(*tvdb_n);i++)
{
	tvdb_tmp2[i]=tvdb_b[i]*tvdb_vg[i];
}
cshift_matrix1D(tvdb_tmp2,(*tvdb_n),1,tvdb_tmp1);
for (i=0;i<(*tvdb_n);i++)
{
	if (vh[i]>0)
	{
		flux1[i]=tvdb_tmp2[i];
	}
	else 
	{
		flux1[i]=tvdb_tmp1[i];
	}
}
//	b1=b-(flux1-cshift(flux1,-1))*dt/2
cshift_matrix1D(flux1,(*tvdb_n),(-1),tvdb_tmp1);
for (i=0;i<(*tvdb_n);i++)
{
	b1[i]=tvdb_b[i]-(flux1[i]-tvdb_tmp1[i])*(*tvdb_dt)/2.0E0;
}
//	do loop start
for (i=0;i<(*tvdb_n);i++)
{
	ip=(i+1)%(*tvdb_n);
	ipp=(ip+1)%(*tvdb_n);
	im=(i+(*tvdb_n)-1)%(*tvdb_n);
	v=vh[i];
	if (v>0)
	{
		w=tvdb_vg[i]*b1[i];
		wp=(tvdb_vg[ip]*b1[ip]-w)/2.0E0;
		wm=(w-tvdb_vg[im]*b1[im])/2.0E0;
	}
	else
	{
		w=tvdb_vg[ip]*b1[ip];
		wp=(w-tvdb_vg[ipp]*b1[ipp])/2.0E0;
		wm=(tvdb_vg[i]*b1[i]-w)/2.0E0;
	}
	dw=0;
	if (wm*wp>0)
	{
		dw=2.0E0*wm*wp/(wm+wp);
	}
	tvdb_flux[i]=(w+dw)*(*tvdb_dt);	
}
//	b=b-(flux-cshift(flux,-1))
cshift_matrix1D(tvdb_flux,(*tvdb_n),(-1),tvdb_tmp1);
for (i=0;i<(*tvdb_n);i++)
{
//	tvdb_flux[i]=vh[i];//tvdb_vg[i];
}
for (i=0;i<(*tvdb_n);i++)
{
	tvdb_b[i]=tvdb_b[i]-(tvdb_flux[i]-tvdb_tmp1[i]);
}
//	do loop end

}
