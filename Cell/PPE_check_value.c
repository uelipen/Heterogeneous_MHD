#include <libspe2.h>
#include <pthread.h>

#include <stdio.h>
#include <stdlib.h>

#include "data_type.h"
#include "parameter.h"

#include "array_definition.h"

void PPE_check_value(data_type_t *PPE_checkV_u, data_type_t *PPE_checkV_b, int *PPE_checkV_nx, int *PPE_checkV_ny, int *PPE_checkV_nz)
{

int i,j,k;

  FILE *check_File_u1;
  FILE *check_File_u2;
  FILE *check_File_u3;
  FILE *check_File_u4;
  FILE *check_File_u5;
  FILE *check_File_b1;
  FILE *check_File_b2;
  FILE *check_File_b3;

//      delete the file
  check_File_u1=fopen("CELL_check_u1.dat","w");
  check_File_u2=fopen("CELL_check_u2.dat","w");
  check_File_u3=fopen("CELL_check_u3.dat","w");
  check_File_u4=fopen("CELL_check_u4.dat","w");
  check_File_u5=fopen("CELL_check_u5.dat","w");
  check_File_b1=fopen("CELL_check_b1.dat","w");
  check_File_b2=fopen("CELL_check_b2.dat","w");
  check_File_b3=fopen("CELL_check_b3.dat","w");

  fclose(check_File_u1);
  fclose(check_File_u2);
  fclose(check_File_u3);
  fclose(check_File_u4);
  fclose(check_File_u5);
  fclose(check_File_b1);
  fclose(check_File_b2);
  fclose(check_File_b3);
//
  check_File_u1=fopen("CELL_check_u1.dat","w+");
  check_File_u2=fopen("CELL_check_u2.dat","w+");
  check_File_u3=fopen("CELL_check_u3.dat","w+");
  check_File_u4=fopen("CELL_check_u4.dat","w+");
  check_File_u5=fopen("CELL_check_u5.dat","w+");
  check_File_b1=fopen("CELL_check_b1.dat","w+");
  check_File_b2=fopen("CELL_check_b2.dat","w+");
  check_File_b3=fopen("CELL_check_b3.dat","w+");
//
printf(" in PPE_check_value\n");
        for (k=0;k<(*PPE_checkV_nz);k++)
        {
                for (j=0;j<(*PPE_checkV_ny);j++)
                {
                        for (i=0;i<(*PPE_checkV_nx);i++)
                        {
                        if (sizeof(data_type_t)==8)
                        {
/*
                                fprintf(check_File_u1,"%i	%i	%i	%17.15E\n",i,j,k,PPE_checkV_u[matrix4D(5,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(1-1),i,j,k)]);
                                fprintf(check_File_u2,"%i       %i      %i      %17.15E\n",i,j,k,PPE_checkV_u[matrix4D(5,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(2-1),i,j,k)]);
                                fprintf(check_File_u3,"%i       %i      %i      %17.15E\n",i,j,k,PPE_checkV_u[matrix4D(5,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(3-1),i,j,k)]);
                                fprintf(check_File_u4,"%i       %i      %i      %17.15E\n",i,j,k,PPE_checkV_u[matrix4D(5,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(4-1),i,j,k)]);
                                fprintf(check_File_u5,"%i       %i      %i      %17.15E\n",i,j,k,PPE_checkV_u[matrix4D(5,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(5-1),i,j,k)]);
                                fprintf(check_File_b1,"%i       %i      %i      %17.15E\n",i,j,k,PPE_checkV_b[matrix4D(3,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(1-1),i,j,k)]);
                                fprintf(check_File_b2,"%i       %i      %i      %17.15E\n",i,j,k,PPE_checkV_b[matrix4D(3,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(2-1),i,j,k)]);
                                fprintf(check_File_b3,"%i       %i      %i      %17.15E\n",i,j,k,PPE_checkV_b[matrix4D(3,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(3-1),i,j,k)]);
*/
//
                                fprintf(check_File_u1,"%17.15E\n",PPE_checkV_u[matrix4D(5,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(1-1),i,j,k)]);
                                fprintf(check_File_u2,"%17.15E\n",PPE_checkV_u[matrix4D(5,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(2-1),i,j,k)]);
                                fprintf(check_File_u3,"%17.15E\n",PPE_checkV_u[matrix4D(5,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(3-1),i,j,k)]);
                                fprintf(check_File_u4,"%17.15E\n",PPE_checkV_u[matrix4D(5,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(4-1),i,j,k)]);
                                fprintf(check_File_u5,"%17.15E\n",PPE_checkV_u[matrix4D(5,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(5-1),i,j,k)]);
                                fprintf(check_File_b1,"%17.15E\n",PPE_checkV_b[matrix4D(3,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(1-1),i,j,k)]);
                                fprintf(check_File_b2,"%17.15E\n",PPE_checkV_b[matrix4D(3,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(2-1),i,j,k)]);
                                fprintf(check_File_b3,"%17.15E\n",PPE_checkV_b[matrix4D(3,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(3-1),i,j,k)]);
//
                        }
                        else if (sizeof(data_type_t)==4)
                        {
                                fprintf(check_File_u1,"%8.6E\n",PPE_checkV_u[matrix4D(5,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(1-1),i,j,k)]);
                                fprintf(check_File_u2,"%8.6E\n",PPE_checkV_u[matrix4D(5,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(2-1),i,j,k)]);
                                fprintf(check_File_u3,"%8.6E\n",PPE_checkV_u[matrix4D(5,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(3-1),i,j,k)]);
                                fprintf(check_File_u4,"%8.6E\n",PPE_checkV_u[matrix4D(5,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(4-1),i,j,k)]);
                                fprintf(check_File_u5,"%8.6E\n",PPE_checkV_u[matrix4D(5,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(5-1),i,j,k)]);
                                fprintf(check_File_b1,"%8.6E\n",PPE_checkV_b[matrix4D(3,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(1-1),i,j,k)]);
                                fprintf(check_File_b2,"%8.6E\n",PPE_checkV_b[matrix4D(3,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(2-1),i,j,k)]);
                                fprintf(check_File_b3,"%8.6E\n",PPE_checkV_b[matrix4D(3,(*PPE_checkV_nx),(*PPE_checkV_ny),(*PPE_checkV_nz),(3-1),i,j,k)]);
                        }
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

return;
}

