#include <stdio.h>
#include <stdlib.h>

#include "array_definition.h"
#include "value_spu.h"

#include "funclist_spu.h"


void tvdb(float *tvdb_flux, float *tvdb_b, float *tvdb_vg, int *tvdb_n, float *tvdb_dt)
{

int i,j,k;

float vh[*tvdb_n];
float flux1[*tvdb_n];
float b1[*tvdb_n];

float tvdb_tmp1[*tvdb_n]; 
float tvdb_tmp2[*tvdb_n];

int ip,ipp,im;
float v,w,wp,wm,dw;

//	vh=(vg+cshift(vg,1))/2
cshift_matrix1D(tvdb_vg,(*tvdb_n),1,tvdb_tmp1);
abXYc_matrix1D(vh,tvdb_vg,tvdb_tmp1,0.5,0.5,(*tvdb_n),0.0);
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
	b1[i]=tvdb_b[i]-(flux1[i]-tvdb_tmp1[i])*(*tvdb_dt)/2.0;
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
		wp=(tvdb_vg[ip]*b1[ip]-w)/2.0;
		wm=(w-tvdb_vg[im]*b1[im])/2.0;
	}
	else
	{
		w=tvdb_vg[ip]*b1[ip];
		wp=(w-tvdb_vg[ipp]*b1[ipp])/2.0;
		wm=(tvdb_vg[i]*b1[i]-w)/2.0;
	}
	dw=0;
	if (wm*wp>0)
	{
		dw=2.0*wm*wp/(wm+wp);
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
