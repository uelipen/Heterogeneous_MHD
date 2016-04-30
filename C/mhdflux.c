#include <math.h>

#include "array_definition.h"

#include "funclist_spu.h"

void find_max(float *max_array, int *max_number, float *max_value);

void mhdflux(float *mhdflux_v, float *mhdflux_c, float *mhdflux_u, float *mhdflux_b, int *mhdflux_n)
{

float vx[*mhdflux_n];
float ps[*mhdflux_n];
float p[*mhdflux_n];
float temp_c[*mhdflux_n];

float gamma;
gamma=5.0/3.0;
int i,ii;

for (i=0;i<(*mhdflux_n);i++)
{
	vx[i]=mhdflux_u[matrix2D(5,(*mhdflux_n),(2-1),i)]/mhdflux_u[matrix2D(5,(*mhdflux_n),(1-1),i)];
}
//
for (i=0;i<(*mhdflux_n);i++)
{
	ps[i]=(mhdflux_u[matrix2D(5,(*mhdflux_n),(5-1),i)]-(mhdflux_u[matrix2D(5,(*mhdflux_n),(2-1),i)]*mhdflux_u[matrix2D(5,(*mhdflux_n),(2-1),i)]+mhdflux_u[matrix2D(5,(*mhdflux_n),(3-1),i)]*mhdflux_u[matrix2D(5,(*mhdflux_n),(3-1),i)]+mhdflux_u[matrix2D(5,(*mhdflux_n),(4-1),i)]*mhdflux_u[matrix2D(5,(*mhdflux_n),(4-1),i)])/mhdflux_u[matrix2D(5,(*mhdflux_n),(1-1),i)]/2.0)*(gamma-1.0)+(2.0-gamma)*(mhdflux_b[matrix2D(3,(*mhdflux_n),(1-1),i)]*mhdflux_b[matrix2D(3,(*mhdflux_n),(1-1),i)]+mhdflux_b[matrix2D(3,(*mhdflux_n),(2-1),i)]*mhdflux_b[matrix2D(3,(*mhdflux_n),(2-1),i)]+mhdflux_b[matrix2D(3,(*mhdflux_n),(3-1),i)]*mhdflux_b[matrix2D(3,(*mhdflux_n),(3-1),i)])/2.0;
}
//
for (i=0;i<(*mhdflux_n);i++)
{
	mhdflux_v[matrix2D(5,(*mhdflux_n),(1-1),i)]=mhdflux_u[matrix2D(5,(*mhdflux_n),(2-1),i)];
	mhdflux_v[matrix2D(5,(*mhdflux_n),(2-1),i)]=mhdflux_u[matrix2D(5,(*mhdflux_n),(2-1),i)]*vx[i]+ps[i]-mhdflux_b[matrix2D(3,(*mhdflux_n),(1-1),i)]*mhdflux_b[matrix2D(3,(*mhdflux_n),(1-1),i)];
	mhdflux_v[matrix2D(5,(*mhdflux_n),(3-1),i)]=mhdflux_u[matrix2D(5,(*mhdflux_n),(3-1),i)]*vx[i]-mhdflux_b[matrix2D(3,(*mhdflux_n),(2-1),i)]*mhdflux_b[matrix2D(3,(*mhdflux_n),(1-1),i)];
	mhdflux_v[matrix2D(5,(*mhdflux_n),(4-1),i)]=mhdflux_u[matrix2D(5,(*mhdflux_n),(4-1),i)]*vx[i]-mhdflux_b[matrix2D(3,(*mhdflux_n),(3-1),i)]*mhdflux_b[matrix2D(3,(*mhdflux_n),(1-1),i)];
	mhdflux_v[matrix2D(5,(*mhdflux_n),(5-1),i)]=(mhdflux_u[matrix2D(5,(*mhdflux_n),(5-1),i)]+ps[i])*vx[i]-mhdflux_b[matrix2D(3,(*mhdflux_n),(1-1),i)]*(mhdflux_b[matrix2D(3,(*mhdflux_n),(1-1),i)]*mhdflux_u[matrix2D(5,(*mhdflux_n),(2-1),i)]+mhdflux_b[matrix2D(3,(*mhdflux_n),(2-1),i)]*mhdflux_u[matrix2D(5,(*mhdflux_n),(3-1),i)]+mhdflux_b[matrix2D(3,(*mhdflux_n),(3-1),i)]*mhdflux_u[matrix2D(5,(*mhdflux_n),(4-1),i)])/mhdflux_u[matrix2D(5,(*mhdflux_n),(1-1),i)];
}
for (i=0;i<(*mhdflux_n);i++)
{
	p[i]=ps[i]-(mhdflux_b[matrix2D(3,(*mhdflux_n),(1-1),i)]*mhdflux_b[matrix2D(3,(*mhdflux_n),(1-1),i)]+mhdflux_b[matrix2D(3,(*mhdflux_n),(2-1),i)]*mhdflux_b[matrix2D(3,(*mhdflux_n),(2-1),i)]+mhdflux_b[matrix2D(3,(*mhdflux_n),(3-1),i)]*mhdflux_b[matrix2D(3,(*mhdflux_n),(3-1),i)])/2.0;
}
for (i=0;i<(*mhdflux_n);i++)
{
	temp_c[i]=fabs(vx[i])+sqrt(fabs((mhdflux_b[matrix2D(3,(*mhdflux_n),(1-1),i)]*mhdflux_b[matrix2D(3,(*mhdflux_n),(1-1),i)]+mhdflux_b[matrix2D(3,(*mhdflux_n),(2-1),i)]*mhdflux_b[matrix2D(3,(*mhdflux_n),(2-1),i)]+mhdflux_b[matrix2D(3,(*mhdflux_n),(3-1),i)]*mhdflux_b[matrix2D(3,(*mhdflux_n),(3-1),i)]+gamma*p[i])/mhdflux_u[matrix2D(5,(*mhdflux_n),(1-1),i)]));
}
find_max(temp_c,mhdflux_n,mhdflux_c);	
if ((*mhdflux_c)>0)
{		
	for (i=0;i<(*mhdflux_n);i++)
	{
		for (ii=0;ii<5;ii++)
		{
			mhdflux_v[matrix2D(5,(*mhdflux_n),ii,i)]=mhdflux_v[matrix2D(5,(*mhdflux_n),ii,i)]/(*mhdflux_c);
		}	
	}
}

}

void find_max(float *max_array, int *max_number, float *max_value)
{
int i;
float temp_max;
temp_max=0.0;
for (i=0; i<(*max_number); i++)
{
        if (max_array[i]>temp_max) temp_max=max_array[i];
}
*max_value=temp_max;
//return;
}
