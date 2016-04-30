#include <stdio.h>
#include <stdlib.h>

#include "array_definition.h"
#include "value_spu.h"

void check_value(float *checkV_u, float *checkV_b, int *checkV_nx, int *checkV_ny, int *checkV_nz, float *checkV_update_u, float *checkV_update_b, int checkV_iter)
{
int i,j,k;

int check_num=50;

  FILE *check_File_u1;
  FILE *check_File_u2;
  FILE *check_File_u3;
  FILE *check_File_u4;
  FILE *check_File_u5;
  FILE *check_File_b1;
  FILE *check_File_b2;
  FILE *check_File_b3;

if (checkV_iter==check_num)
{
//      delete the file
  check_File_u1=fopen("C_version_check_u1.dat","w");
  check_File_u2=fopen("C_version_check_u2.dat","w");
  check_File_u3=fopen("C_version_check_u3.dat","w");
  check_File_u4=fopen("C_version_check_u4.dat","w");
  check_File_u5=fopen("C_version_check_u5.dat","w");
  check_File_b1=fopen("C_version_check_b1.dat","w");
  check_File_b2=fopen("C_version_check_b2.dat","w");
  check_File_b3=fopen("C_version_check_b3.dat","w");

  fclose(check_File_u1);
  fclose(check_File_u2);
  fclose(check_File_u3);
  fclose(check_File_u4);
  fclose(check_File_u5);
  fclose(check_File_b1);
  fclose(check_File_b2);
  fclose(check_File_b3);
//
  printf("%i\n",checkV_iter);
  check_File_u1=fopen("C_version_check_u1.dat","w+");
  check_File_u2=fopen("C_version_check_u2.dat","w+");
  check_File_u3=fopen("C_version_check_u3.dat","w+");
  check_File_u4=fopen("C_version_check_u4.dat","w+");
  check_File_u5=fopen("C_version_check_u5.dat","w+");
  check_File_b1=fopen("C_version_check_b1.dat","w+");
  check_File_b2=fopen("C_version_check_b2.dat","w+");
  check_File_b3=fopen("C_version_check_b3.dat","w+");

	for (k=0;k<(*checkV_nz);k++)
	{
		for (j=0;j<(*checkV_ny);j++)
		{
			for (i=0;i<(*checkV_nx);i++)
			{
				fprintf(check_File_u1,"%8.6e\n",checkV_u[matrix4D(5,(*checkV_nx),(*checkV_ny),(*checkV_nz),(1-1),i,j,k)]);
				fprintf(check_File_u2,"%8.6e\n",checkV_u[matrix4D(5,(*checkV_nx),(*checkV_ny),(*checkV_nz),(2-1),i,j,k)]);
				fprintf(check_File_u3,"%8.6e\n",checkV_u[matrix4D(5,(*checkV_nx),(*checkV_ny),(*checkV_nz),(3-1),i,j,k)]);
				fprintf(check_File_u4,"%8.6e\n",checkV_u[matrix4D(5,(*checkV_nx),(*checkV_ny),(*checkV_nz),(4-1),i,j,k)]);
				fprintf(check_File_u5,"%8.6e\n",checkV_u[matrix4D(5,(*checkV_nx),(*checkV_ny),(*checkV_nz),(5-1),i,j,k)]);
                                fprintf(check_File_b1,"%8.6e\n",checkV_b[matrix4D(3,(*checkV_nx),(*checkV_ny),(*checkV_nz),(1-1),i,j,k)]);
                                fprintf(check_File_b2,"%8.6e\n",checkV_b[matrix4D(3,(*checkV_nx),(*checkV_ny),(*checkV_nz),(2-1),i,j,k)]);
                                fprintf(check_File_b3,"%8.6e\n",checkV_b[matrix4D(3,(*checkV_nx),(*checkV_ny),(*checkV_nz),(3-1),i,j,k)]);
			}
		}
	}
  fclose(check_File_u1);
  fclose(check_File_u2);
  fclose(check_File_u3);
  fclose(check_File_u4);
  fclose(check_File_u5);
  fclose(check_File_b1);
  fclose(check_File_b2);
  fclose(check_File_b3);

}
}
