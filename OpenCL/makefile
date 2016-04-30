
#machine_selection=ant
machine_selection=muskrat
#machine_selection=gamma 		# it doesn't work yet
#
#
exec = mhd 
#
CC = g++
CFLAGS = 
LDLIBS = -lOpenCL -lSDKUtil

ifeq ($(machine_selection),ant)
LDPATH =  -L/usr/ati-stream-sdk-v2.0-beta4-lnx64/lib/x86_64
LDINC = -I/usr/ati-stream-sdk-v2.0-beta4-lnx64/include 
endif

ifeq ($(machine_selection),muskrat)
LDPATH =  -L/usr/lib64/ -L/cita/h/home-1/bpang/NVIDIA_GPU_Computing_SDK/lib
LDINC = -I/cita/h/home-1/bpang/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc 
endif

ifeq ($(machine_selection),gamma)
LDPATH =  -L/cita/h/home-1/bpang/gamma_opencl_lib
LDINC = -I/cita/h/home-1/bpang/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc
endif

objects = mhd.o

$(exec): $(objects)
	$(CC) $(CFLAGS) $(objects) -o $@ $(LDPATH) $(LDLIBS) $(LDINC)

mhd.o: mhd.cpp 
	$(CC) $(CFLAGS) -c $< -o $@ $(LDPATH) $(LDLIBS) $(LDINC)

mhd.o: mhd.hpp

.PHONY: clean
clean:
	rm -f $(exec) *.o

backup:
	tar -cvf backup_OpenCL_mhd.tar *.cpp *.hpp *.cl makefile ; cp backup_OpenCL_mhd.tar ~/opencl_MHD/BACKUP/

