


DIRS		:= spu

PROGRAM_ppu	:= mhd_ppu 

IMPORTS		= spu/lib_mhd_spu.a -lspe2 -lpthread -lm

CC_OPT_LEVEL    = -O0 -g

ifdef CELL_TOP
        include $(CELL_TOP)/buildutils/make.footer
else
        include ../../../../buildutils/make.footer
endif

backup:
	rm backup_multiCore_cell_mhd.tar; tar -cvf backup_multiCore_cell_mhd.tar *.c *.h Makefile spu/*.c spu/*.h spu/Makefile; cp backup_multiCore_cell_mhd.tar ~/BACKUP/cell/
