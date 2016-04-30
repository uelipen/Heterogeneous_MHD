#ifndef INCLUDEFUNC_SPU_H
#define INCLUDEFUNC_SPU_H

void multiple_itself_matrix4D(float *mim4_ub, int mim4_nub, int mim4_nx, int mim4_ny, int mim4_nz, float mim4_multiple)
{

int ii,i,j,k;
for (k=0;k<mim4_nz;k++)
{
	for (j=0;j<mim4_ny;j++)
	{
		for (i=0;i<mim4_nx;i++)
		{
			for (ii=0;ii<mim4_nub;ii++)
			{
				mim4_ub[matrix4D(mim4_nub,mim4_nx,mim4_ny,mim4_nz,ii,i,j,k)]=mim4_multiple*mim4_ub[matrix4D(mim4_nub,mim4_nx,mim4_ny,mim4_nz,ii,i,j,k)];
			}
		}
	}
}

return;
}

void multiple_itself_matrix2D(float *mim2_ub, int mim2_nub, int mim2_nx, float mim2_multiple)
{

int ii,i;
for (i=0;i<mim2_nx;i++)
{
	for (ii=0;ii<mim2_nub;ii++)
        {
		mim2_ub[matrix2D(mim2_nub,mim2_nx,ii,i)]=mim2_multiple*mim2_ub[matrix2D(mim2_nub,mim2_nx,ii,i)];
     	}
}
return;
}

void cshift_matrix1D(float *cm1_ub, int cm1_nx, int cm1_logical, float *cm1_out_ub)
{
int i;
int imp,imm;

for (i=0;i<(cm1_nx);i++)
{
	imm=(i+(cm1_nx)-1)%(cm1_nx);
	imp=(i+1)%(cm1_nx);
	if (cm1_logical==1)
	{
		cm1_out_ub[i]=cm1_ub[imp];

	}
	else if (cm1_logical==-1)
	{
                cm1_out_ub[i]=cm1_ub[imm];
	}
}

return;
}

void extract_matrix2Dto1D(float *em2to1_ub, int em2to1_nub, int em2to1_nx, int em2to1_ii, float *em2to1_out_ub)
{
int i;

for (i=0;i<em2to1_nx;i++)
{
	em2to1_out_ub[i]=em2to1_ub[matrix2D(em2to1_nub,em2to1_nx,em2to1_ii,i)];
}
return;
}

void abXYc_matrix2D(float *abXYcm2D_L, float *abXYcm2D_Ra, float *abXYcm2D_Rb, float abXYcm2D_a, float abXYcm2D_b, int abXYcm2D_nub, int abXYcm2D_nx, float abXYcm2D_c)
{
int i,ii;
for (i=0;i<abXYcm2D_nx;i++)
{
	for (ii=0;ii<abXYcm2D_nub;ii++)
	{
		abXYcm2D_L[matrix2D(abXYcm2D_nub,abXYcm2D_nx,ii,i)]=abXYcm2D_a*abXYcm2D_Ra[matrix2D(abXYcm2D_nub,abXYcm2D_nx,ii,i)]+abXYcm2D_b*abXYcm2D_Rb[matrix2D(abXYcm2D_nub,abXYcm2D_nx,ii,i)]+abXYcm2D_c;
	}
}
return;
}

void abXYc_matrix1D(float *abXYcm1D_L, float *abXYcm1D_Ra, float *abXYcm1D_Rb, float abXYcm1D_a, float abXYcm1D_b, int abXYcm1D_nx, float abXYcm1D_c)
{
int i;
for (i=0;i<abXYcm1D_nx;i++)
{
	abXYcm1D_L[i]=abXYcm1D_a*abXYcm1D_Ra[i]+abXYcm1D_b*abXYcm1D_Rb[i]+abXYcm1D_c;
}
return;
}

void cshift_matrix2D_dim2(float *cm2d2_ub, int cm2d2_nub, int cm2d2_nx, int cm2d2_logical, float *cm2d2_out_ub)
{
int i,ii;
int imm,imp;
for (i=0;i<cm2d2_nx;i++)
{
        imm=(i+(cm2d2_nx)-1)%(cm2d2_nx);
        imp=(i+1)%(cm2d2_nx);
	for (ii=0;ii<cm2d2_nub;ii++)
	{
		if (cm2d2_logical==1)
		{
			cm2d2_out_ub[matrix2D(cm2d2_nub,cm2d2_nx,ii,i)]=cm2d2_ub[matrix2D(cm2d2_nub,cm2d2_nx,ii,imp)];
		}
		else if (cm2d2_logical==-1)
		{
                        cm2d2_out_ub[matrix2D(cm2d2_nub,cm2d2_nx,ii,i)]=cm2d2_ub[matrix2D(cm2d2_nub,cm2d2_nx,ii,imm)];
		}
	}
}
return;
}

void limiter_matrix2D(float *lm2D_ubA, float *lm2D_ubB, int lm2D_nub, int lm2D_nx, float *lm2D_out_ub)
{
int i,ii;
for (i=0;i<lm2D_nx;i++)
{
	for (ii=0;ii<lm2D_nub;ii++)
	{
		if (lm2D_ubA[matrix2D(lm2D_nub,lm2D_nx,ii,i)]*lm2D_ubB[matrix2D(lm2D_nub,lm2D_nx,ii,i)]>0)
		{
			lm2D_out_ub[matrix2D(lm2D_nub,lm2D_nx,ii,i)]=2.0*lm2D_ubA[matrix2D(lm2D_nub,lm2D_nx,ii,i)]*lm2D_ubB[matrix2D(lm2D_nub,lm2D_nx,ii,i)]/(lm2D_ubA[matrix2D(lm2D_nub,lm2D_nx,ii,i)]+lm2D_ubB[matrix2D(lm2D_nub,lm2D_nx,ii,i)]);
		}
	}
}

return;
}

#endif
