exec = mhd_exec

FC = gcc

FFLAGS = 

LDLIBS = -lm


objects = check_value.o	\
	check_update.o	\
	calcfl.o	\
	fluidx.o mhdflux.o mhd.o copy_matrix.o tvd1.o tvdb.o init.o \
	advectbyzx.o	\
	advectbyzxA1.o	\
	advectbyzxA2.o	\
	advectbyzxB1.o	\
	advectbyzxB2.o	\
	transpose12.o	\
	transpose13.o	

$(exec): $(objects)
	$(FC) $(FFLAGS) $(objects) -o $@ $(LDPATH) $(LDLIBS)

calcfl.o: calcfl.c
	$(FC) $(FFLAGS) -c $< -o $@

check_update.o: check_update.c
	$(FC) $(FFLAGS) -c $< -o $@

check_value.o: check_value.c
	$(FC) $(FFLAGS) -c $< -o $@

transpose13.o: transpose13.c
	$(FC) $(FFLAGS) -c $< -o $@

transpose12.o: transpose12.c
	$(FC) $(FFLAGS) -c $< -o $@

advectbyzxB2.o: advectbyzxB2.c
	$(FC) $(FFLAGS) -c $< -o $@

advectbyzxB1.o: advectbyzxB1.c
	$(FC) $(FFLAGS) -c $< -o $@

advectbyzxA2.o: advectbyzxA2.c
	$(FC) $(FFLAGS) -c $< -o $@

advectbyzxA1.o: advectbyzxA1.c
	$(FC) $(FFLAGS) -c $< -o $@

advectbyzx.o: advectbyzx.c
	$(FC) $(FFLAGS) -c $< -o $@

fluidx.o: fluidx.c
	$(FC) $(FFLAGS) -c $< -o $@

mhdflux.o: mhdflux.c
	$(FC) $(FFLAGS) -c $< -o $@

mhd.o: mhd.c
	$(FC) $(FFLAGS) -c $< -o $@

copy_matrix.o: copy_matrix.c
	$(FC) $(FFLAGS) -c $< -o $@

tvdb.o: tvdb.c
	$(FC) $(FFLAGS) -c $< -o $@

tvd1.o: tvd1.c
	$(FC) $(FFLAGS) -c $< -o $@

init.o: init.c
	$(FC) $(FFLAGS) -c $< -o $@


.PHONY: clean
clean:
	rm -f $(exec) *.o

backup:
	tar -cvf backup_C_mhd.tar *.c *.h makefile ; cp backup_C_mhd.tar ~/BACKUP/

