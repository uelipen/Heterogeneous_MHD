#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#include <stdio.h>
#include <stdlib.h>

#include "../data_type.h"
#include "../parameter.h"

#include "array_definition.h"
#include "value_spu.h"

void check_value(data_type_t *checkV_u, data_type_t *checkV_b, int *checkV_nx, int *checkV_ny, int *checkV_nz, int *checkV_Y_location, int *checkV_Z_location, int checkV_iter, int checkV_SPE_id, int checkV_tag_id)
{

volatile data_type_t spu_checkV_u[5*(*checkV_nx)];
volatile data_type_t spu_checkV_b[3*(*checkV_nx)];

int i,j,k;
int j_global,k_global;

int Y_inter=box_ny/num_SPE_Y;
int Z_inter=box_nz/num_SPE_Z;

int CELL_PER_BLOCK;
CELL_PER_BLOCK = (*checkV_nx);
int u_memSize_PER_BLOCK;
u_memSize_PER_BLOCK = 5* CELL_PER_BLOCK;
int b_memSize_PER_BLOCK;
b_memSize_PER_BLOCK = 3* CELL_PER_BLOCK;

if (checkV_SPE_id==0)
{
int check_num=10;

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
//
  char name_file_id[5];
  char file_name_u1[50];
  char file_name_u2[50];
  char file_name_u3[50];
  char file_name_u4[50];
  char file_name_u5[50];
  char file_name_b1[50];
  char file_name_b2[50];
  char file_name_b3[50];
//
  printf("%i\n",checkV_iter);
  int n_file;
  n_file=sprintf (name_file_id, "00%d", checkV_SPE_id);
  sprintf(file_name_u1,"CELL_version_check_u1_%s.dat",name_file_id);
  sprintf(file_name_u2,"CELL_version_check_u2_%s.dat",name_file_id);
  sprintf(file_name_u3,"CELL_version_check_u3_%s.dat",name_file_id);
  sprintf(file_name_u4,"CELL_version_check_u4_%s.dat",name_file_id);
  sprintf(file_name_u5,"CELL_version_check_u5_%s.dat",name_file_id);
  sprintf(file_name_b1,"CELL_version_check_b1_%s.dat",name_file_id);
  sprintf(file_name_b2,"CELL_version_check_b2_%s.dat",name_file_id);
  sprintf(file_name_b3,"CELL_version_check_b3_%s.dat",name_file_id);

printf("id is %i\n",checkV_SPE_id);
/*
printf("u1 file is %s\n",file_name_u1);
printf("u2 file is %s\n",file_name_u2);
printf("u3 file is %s\n",file_name_u3);
printf("u4 file is %s\n",file_name_u4);
printf("u5 file is %s\n",file_name_u5);
printf("b1 file is %s\n",file_name_b1);
printf("b2 file is %s\n",file_name_b2);
printf("b3 file is %s\n",file_name_b3);
*/
  barrier(checkV_SPE_id);
//
  check_File_u1=fopen(file_name_u1,"w+");
  check_File_u2=fopen(file_name_u2,"w+");
  check_File_u3=fopen(file_name_u3,"w+");
  check_File_u4=fopen(file_name_u4,"w+");
  check_File_u5=fopen(file_name_u5,"w+");
  check_File_b1=fopen(file_name_b1,"w+");
  check_File_b2=fopen(file_name_b2,"w+");
  check_File_b3=fopen(file_name_b3,"w+");
  barrier(checkV_SPE_id);
printf(" in check_value\n");
//      delete the file
/*
  fclose(check_File_u1);
  fclose(check_File_u2);
  fclose(check_File_u3);
  fclose(check_File_u4);
  fclose(check_File_u5);
  fclose(check_File_b1);
  fclose(check_File_b2);
  fclose(check_File_b3);
//
  check_File_u1=fopen(file_name_u1,"w+");
  check_File_u2=fopen(file_name_u2,"w+");
  check_File_u3=fopen(file_name_u3,"w+");
  check_File_u4=fopen(file_name_u4,"w+");
  check_File_u5=fopen(file_name_u5,"w+");
  check_File_b1=fopen(file_name_b1,"w+");
  check_File_b2=fopen(file_name_b2,"w+");
  check_File_b3=fopen(file_name_b3,"w+");
*/
//

  for (k=0;k<(*checkV_nz);k++)
  {
        for (j=0;j<(*checkV_ny);j++)
        {
//	assign global index
		k_global=k+(*checkV_Z_location)*Z_inter;
		j_global=j+(*checkV_Y_location)*Y_inter;
//	get u
                spu_mfcdma32((void *)(spu_checkV_u),
		(unsigned int)(checkV_u+matrix4D(5,(box_nx),(box_ny),(box_nz),(1-1),(1-1),(j_global),(k_global))), 
		u_memSize_PER_BLOCK * sizeof(data_type_t), checkV_tag_id, MFC_GETB_CMD);
//	get b
                spu_mfcdma32((void *)(spu_checkV_b),
                (unsigned int)(checkV_b+matrix4D(3,(box_nx),(box_ny),(box_nz),(1-1),(1-1),(j_global),(k_global))),
                b_memSize_PER_BLOCK * sizeof(data_type_t), checkV_tag_id, MFC_GET_CMD);
//	Wait for the DMA to complete
                (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);
//	now you can print 
//
                for (i=0;i<(*checkV_nx);i++)
                {
                        if (sizeof(data_type_t)==4)
                        {
                        fprintf(check_File_u1,"%8.6E\n",spu_checkV_u[matrix2D(5,(*checkV_nx),(1-1),i)]);
                        fprintf(check_File_u2,"%8.6E\n",spu_checkV_u[matrix2D(5,(*checkV_nx),(2-1),i)]);
                        fprintf(check_File_u3,"%8.6E\n",spu_checkV_u[matrix2D(5,(*checkV_nx),(3-1),i)]);
                        fprintf(check_File_u4,"%8.6E\n",spu_checkV_u[matrix2D(5,(*checkV_nx),(4-1),i)]);
                        fprintf(check_File_u5,"%8.6E\n",spu_checkV_u[matrix2D(5,(*checkV_nx),(5-1),i)]);
                        fprintf(check_File_b1,"%8.6E\n",spu_checkV_b[matrix2D(3,(*checkV_nx),(1-1),i)]);
                        fprintf(check_File_b2,"%8.6E\n",spu_checkV_b[matrix2D(3,(*checkV_nx),(2-1),i)]);
                        fprintf(check_File_b3,"%8.6E\n",spu_checkV_b[matrix2D(3,(*checkV_nx),(3-1),i)]);
                        }
                        else if (sizeof(data_type_t)==8)
                        {
                        fprintf(check_File_u1,"%17.15E\n",spu_checkV_u[matrix2D(5,(*checkV_nx),(1-1),i)]);
                        fprintf(check_File_u2,"%17.15E\n",spu_checkV_u[matrix2D(5,(*checkV_nx),(2-1),i)]);
                        fprintf(check_File_u3,"%17.15E\n",spu_checkV_u[matrix2D(5,(*checkV_nx),(3-1),i)]);
                        fprintf(check_File_u4,"%17.15E\n",spu_checkV_u[matrix2D(5,(*checkV_nx),(4-1),i)]);
                        fprintf(check_File_u5,"%17.15E\n",spu_checkV_u[matrix2D(5,(*checkV_nx),(5-1),i)]);
                        fprintf(check_File_b1,"%17.15E\n",spu_checkV_b[matrix2D(3,(*checkV_nx),(1-1),i)]);
                        fprintf(check_File_b2,"%17.15E\n",spu_checkV_b[matrix2D(3,(*checkV_nx),(2-1),i)]);
                        fprintf(check_File_b3,"%17.15E\n",spu_checkV_b[matrix2D(3,(*checkV_nx),(3-1),i)]);
printf("fk spu_checkV_u[matrix2D(5,(*checkV_nx),(1-1),i)]: %e\n",spu_checkV_u[matrix2D(5,(*checkV_nx),(1-1),i)]);
                        }
                }

        }
  }

//      Wait for final DMAs to complete before terminating SPU thread.
  (void)spu_mfcstat(MFC_TAG_UPDATE_ALL);

}
}
return;
}

