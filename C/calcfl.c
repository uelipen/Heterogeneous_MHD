#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "array_definition.h"
#include "value_spu.h"

#include "funclist_spu.h"

float max_3num(float *m3_1, float *m3_2, float *m3_3);
float max_2num(float *m2_1, float *m2_2);

float calcfl(float *calcfl_u, float *calcfl_b, int *calcfl_nx, int *calcfl_ny, int *calcfl_nz)
{
int i,j,k;

float gamma;
gamma=5.0/3.0;
float c;
c=0.0;

int ip,jp,kp;
float bx,by,bz,v,ps,p;
float temp1,temp2,temp3;

for (k=0;k<(*calcfl_nz);k++)
{
	for (j=0;j<(*calcfl_ny);j++)
	{
		for (i=0;i<(*calcfl_nx);i++)
		{
//	kp=mod(k,nz)+1
			kp=(k+1)%(*calcfl_nz);	
//	jp=mod(j,ny)+1
			jp=(j+1)%(*calcfl_ny);	
//	ip=mod(i,nx)+1
			ip=(i+1)%(*calcfl_nx);	
//	bx=(b(1,i,j,k)+b(1,ip,j,k))/2.0d0
			bx=(calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(1-1),i,j,k)]+calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(1-1),ip,j,k)])/2.0;
//	by=(b(2,i,j,k)+b(2,i,jp,k))/2.0d0
			by=(calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(2-1),i,j,k)]+calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(2-1),i,jp,k)])/2.0;
//	bz=(b(3,i,j,k)+b(3,i,j,kp))/2.0d0
			bz=(calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(3-1),i,j,k)]+calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(3-1),i,j,kp)])/2.0;
//	v=maxval(abs(u(2:4,i,j,k)/u(1,i,j,k)))
			temp1=fabs(calcfl_u[matrix4D(5,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(2-1),i,j,k)]/calcfl_u[matrix4D(5,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(1-1),i,j,k)]);
			temp2=fabs(calcfl_u[matrix4D(5,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(3-1),i,j,k)]/calcfl_u[matrix4D(5,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(1-1),i,j,k)]);
			temp3=fabs(calcfl_u[matrix4D(5,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(4-1),i,j,k)]/calcfl_u[matrix4D(5,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(1-1),i,j,k)]);
			v=max_3num(&temp1,&temp2,&temp3);
//	ps=(u(5,i,j,k)-sum(u(2:4,i,j,k)**2,1)/u(1,i,j,k)/2.0d0)*(gamma-1.0d0)+(2.0d0-gamma)*sum(b(:,i,j,k)**2,1)/2.0d0
			ps=(calcfl_u[matrix4D(5,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(5-1),i,j,k)]-(calcfl_u[matrix4D(5,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(2-1),i,j,k)]*calcfl_u[matrix4D(5,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(2-1),i,j,k)]+calcfl_u[matrix4D(5,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(3-1),i,j,k)]*calcfl_u[matrix4D(5,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(3-1),i,j,k)]+calcfl_u[matrix4D(5,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(4-1),i,j,k)]*calcfl_u[matrix4D(5,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(4-1),i,j,k)])/calcfl_u[matrix4D(5,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(1-1),i,j,k)]/2.0)*(gamma-1.0)+(2.0-gamma)*(calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(1-1),i,j,k)]*calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(1-1),i,j,k)]+calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(2-1),i,j,k)]*calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(2-1),i,j,k)]+calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(3-1),i,j,k)]*calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(3-1),i,j,k)])/2.0;
//	p=ps-sum(b(:,i,j,k)**2)/2.0d0
			p=ps-(calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(1-1),i,j,k)]*calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(1-1),i,j,k)]+calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(2-1),i,j,k)]*calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(2-1),i,j,k)]+calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(3-1),i,j,k)]*calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(3-1),i,j,k)])/2.0;
//	c=max(c,v+sqrt(abs(  (sum(b(:,i,j,k)**2)*2.0d0+gamma*p)/u(1,i,j,k))))			
			temp1=v+sqrt(fabs((calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(1-1),i,j,k)]*calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(1-1),i,j,k)]+calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(2-1),i,j,k)]*calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(2-1),i,j,k)]+calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(3-1),i,j,k)]*calcfl_b[matrix4D(3,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(3-1),i,j,k)])*2.0+gamma*p)/calcfl_u[matrix4D(5,(*calcfl_nx),(*calcfl_ny),(*calcfl_nz),(1-1),i,j,k)]);
			temp2=max_2num(&c,&temp1);
			c=temp2;
		}
	}
}

temp1=1/c;
return temp1;
}

float max_3num(float *m3_1, float *m3_2, float *m3_3)
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

float max_2num(float *m2_1, float *m2_2)
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
