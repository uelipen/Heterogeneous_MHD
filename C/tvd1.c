#include "array_definition.h"
#include <assert.h>

#include "funclist_spu.h"

void mhdflux(float *mhdflux_v, float *mhdflux_c, float *mhdflux_u, float *mhdflux_b, int *mhdflux_n);

void tvd1(float *tvd1_u, float *tvd1_b, int *tvd1_n, float *tvd1_dt)
{

float wr[5*(*tvd1_n)]; 
float wl[5*(*tvd1_n)]; 
float v[5*(*tvd1_n)];
float fr[5*(*tvd1_n)];
float fl[5*(*tvd1_n)];
float flux[5*(*tvd1_n)];
float u1[5*(*tvd1_n)];

float dfrp[5*(*tvd1_n)];
float dfrm[5*(*tvd1_n)];
float dfr[5*(*tvd1_n)];
float dflp[5*(*tvd1_n)];
float dflm[5*(*tvd1_n)];
float dfl[5*(*tvd1_n)];

float c;

float tmp1[5*(*tvd1_n)];
float tmp2[5*(*tvd1_n)];

int i,j,k,ii;

//	mhdflux
mhdflux(v,&c,tvd1_u,tvd1_b,tvd1_n);
//	wr=u+v
abXYc_matrix2D(wr,tvd1_u,v,1.0,1.0,5,(*tvd1_n),0.0);
//	wl=u-v
abXYc_matrix2D(wl,tvd1_u,v,1.0,-1.0,5,(*tvd1_n),0.0);
//	fr=c*wr
abXYc_matrix2D(fr,wr,wr,c,0.0,5,(*tvd1_n),0.0);
//	fl=cshift(c*wl,1,2)
cshift_matrix2D_dim2(wl,5,(*tvd1_n),1,tmp1);
abXYc_matrix2D(fl,tmp1,tmp1,c,0.0,5,(*tvd1_n),0.0);
//	flux=(fr-fl)/2
abXYc_matrix2D(flux,fr,fl,0.5,-0.5,5,(*tvd1_n),0.0);
//	u1=u-(flux-cshift(flux,-1,2))*dt/2
cshift_matrix2D_dim2(flux,5,(*tvd1_n),-1,tmp1);
abXYc_matrix2D(tmp2,flux,tmp1,1.0,-1.0,5,(*tvd1_n),0.0);
abXYc_matrix2D(u1,tvd1_u,tmp2,1.0,(-1.0*(*tvd1_dt)/2.0),5,(*tvd1_n),0.0);

//	mhdflux
mhdflux(v,&c,u1,tvd1_b,tvd1_n);
//	wr=u1+v
abXYc_matrix2D(wr,u1,v,1.0,1.0,5,(*tvd1_n),0.0);
//	wl=u1-v
abXYc_matrix2D(wl,u1,v,1.0,-1.0,5,(*tvd1_n),0.0);
//	fr=c*wr
abXYc_matrix2D(fr,wr,wr,c,0.0,5,(*tvd1_n),0.0);
//	dfrp=(cshift(fr,1,2)-fr)/2
cshift_matrix2D_dim2(fr,5,(*tvd1_n),1,tmp1);
abXYc_matrix2D(dfrp,tmp1,fr,0.5,-0.5,5,(*tvd1_n),0.0);
//	dfrm=(fr-cshift(fr,-1,2))/2
cshift_matrix2D_dim2(fr,5,(*tvd1_n),-1,tmp2);
abXYc_matrix2D(dfrm,fr,tmp2,0.5,-0.5,5,(*tvd1_n),0.0);
//	dfr=0
//abXYc_matrix2D(dfr,dfr,dfr,0.0,0.0,5,(*tvd1_n),0);
//	where(dfrp*dfrm>0)
//	dfr=2*dfrp*dfrm/(dfrp+dfrm)
//	end where
//limiter_matrix2D(dfrp,dfrm,5,(*tvd1_n),dfr);
for (i=0;i<(*tvd1_n);i++)
{
	for (ii=0;ii<5;ii++)
	{
		dfr[matrix2D(5,(*tvd1_n),ii,i)]=0.0;
		if (dfrp[matrix2D(5,(*tvd1_n),ii,i)]*dfrm[matrix2D(5,(*tvd1_n),ii,i)]>0)
		{
			dfr[matrix2D(5,(*tvd1_n),ii,i)]=2.0*dfrp[matrix2D(5,(*tvd1_n),ii,i)]*dfrm[matrix2D(5,(*tvd1_n),ii,i)]/(dfrp[matrix2D(5,(*tvd1_n),ii,i)]+dfrm[matrix2D(5,(*tvd1_n),ii,i)]);
		}
	}
}
//	fl=cshift(c*wl,1,2)
cshift_matrix2D_dim2(wl,5,(*tvd1_n),1,tmp2);
abXYc_matrix2D(fl,tmp2,tmp2,c,0.0,5,(*tvd1_n),0);
//	dflp=(fl-cshift(fl,1,2))/2
cshift_matrix2D_dim2(fl,5,(*tvd1_n),1,tmp1);
abXYc_matrix2D(dflp,fl,tmp1,0.5,-0.5,5,(*tvd1_n),0.0);
//	dflm=(cshift(fl,-1,2)-fl)/2
cshift_matrix2D_dim2(fl,5,(*tvd1_n),-1,tmp2);
abXYc_matrix2D(dflm,tmp2,fl,0.5,-0.5,5,(*tvd1_n),0.0);
//	dfl=0
//abXYc_matrix2D(dfl,dfl,dfl,0.0,0.0,5,(*tvd1_n),0.0);
//	where(dflp*dflm>0)
//	dfl=2*dflp*dflm/(dflp+dflm)
//	end where
//limiter_matrix2D(dflp,dflm,5,(*tvd1_n),dfl);
for (i=0;i<(*tvd1_n);i++)
{
        for (ii=0;ii<5;ii++)
        {
                dfl[matrix2D(5,(*tvd1_n),ii,i)]=0.0;
                if (dflp[matrix2D(5,(*tvd1_n),ii,i)]*dflm[matrix2D(5,(*tvd1_n),ii,i)]>0)
                {
                        dfl[matrix2D(5,(*tvd1_n),ii,i)]=2.0*dflp[matrix2D(5,(*tvd1_n),ii,i)]*dflm[matrix2D(5,(*tvd1_n),ii,i)]/(dflp[matrix2D(5,(*tvd1_n),ii,i)]+dflm[matrix2D(5,(*tvd1_n),ii,i)]);
                }
        }
}

//	flux=(fr-fl+(dfr-dfl))/2
abXYc_matrix2D(tmp1,fr,fl,1.0,-1.0,5,(*tvd1_n),0.0);
abXYc_matrix2D(tmp2,dfr,dfl,1.0,-1.0,5,(*tvd1_n),0.0);
abXYc_matrix2D(flux,tmp1,tmp2,0.5,0.5,5,(*tvd1_n),0.0);
//	u=u-(flux-cshift(flux,-1,2))*dt
cshift_matrix2D_dim2(flux,5,(*tvd1_n),-1,tmp1);
abXYc_matrix2D(tmp2,flux,tmp1,1.0,-1.0,5,(*tvd1_n),0.0);
abXYc_matrix2D(tvd1_u,tvd1_u,tmp2,1.0,(-1.0*(*tvd1_dt)),5,(*tvd1_n),0.0);
//	finished tvd1
}
