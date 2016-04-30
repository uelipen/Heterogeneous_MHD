exec = Cudamhd

# Define Fortran compiler
CC = gcc
LDPATH = -L /cita/h/home-1/bpang/CUDA/cuda/lib
LDLIBS = -lcudart #-lcrt
FFLAGS = #-r8 #-g -CB -fpe0
#
cuda_C = nvcc
cuda_FFLAGS = #-deviceemu#-g# -O3
cuda_LDLIBS = #-lc
#FC= g95

objects = mhd.o	\
	 init.o	\
	 cuda_main.o	\
	 cuda_cfl.o	\
	 cuda_fluidx.o	\
	 cuda_advectbyzx.o	\
	 cuda_transpose12.o	\
	 cuda_transpose13.o

$(exec): $(objects)
	$(CC) $(FFLAGS) $(objects) -o $@ $(LDPATH) $(LDLIBS) 

mhd.o: mhd.c
	$(CC) $(FFLAGS) -c $< -o $@

init.o: init.c
	$(CC) $(FFLAGS) -c $< -o $@

cuda_main.o: cuda_main.cu
	$(cuda_C) $(cuda_FFLAGS) -c $< -o $@ $(cuda_LDLIBS)

cuda_cfl.o: cuda_cfl.cu
	$(cuda_C) $(cuda_FFLAGS) -c $< -o $@ $(cuda_LDLIBS)

cuda_fluidx.o: cuda_fluidx.cu
	$(cuda_C) $(cuda_FFLAGS) -c $< -o $@ $(cuda_LDLIBS)

cuda_advectbyzx.o: cuda_advectbyzx.cu
	$(cuda_C) $(cuda_FFLAGS) -c $< -o $@ $(cuda_LDLIBS)

cuda_transpose12.o: cuda_transpose12.cu
	$(cuda_C) $(cuda_FFLAGS) -c $< -o $@ $(cuda_LDLIBS)

cuda_transpose13.o: cuda_transpose13.cu
	$(cuda_C) $(cuda_FFLAGS) -c $< -o $@ $(cuda_LDLIBS)

.PHONY: clean
clean: 
	rm -f $(exec) *.o *.linkinfo

backup:
	tar -cvf backup_cuda_mhd.tar *.c *.cu *.h makefile ; cp backup_cuda_mhd.tar ~/cuda_MHD/BACKUP/
# tar -xvf
