###################################################################
#
# freeman.justin@gmail.com 
#
##################################################################

CC=	gcc

CSRC=	./src/

CFLAGS=	-O3 -g -Wall 
#CFLAGS=	-O3 -g -fPIC -Wall

NNLIB=	./lib/nn

INC=	-I./include \
	-I$(NNLIB)


#LFLAGS= -lnetcdf -lxml2 
LFLAGS= -ludunits2 -lnetcdf $(NNLIB)/libnn.a

COBJ=	$(CSRC)main.o \
	$(CSRC)jutil.o \
	$(CSRC)fail.o \
	$(CSRC)spheriq_dist.o \
	$(CSRC)malloc_arrays.o \
	$(CSRC)netcdfIO.o \
	$(CSRC)interp.o \
	$(CSRC)calcalcs.o \
	$(CSRC)utCalendar2_cal.o \
	$(CSRC)process_roms.o \
	$(CSRC)process_tides.o \
	$(CSRC)process_waves.o \
	$(CSRC)process_reference_levels.o \
	$(CSRC)add.o \
	$(CSRC)cli.o \
	$(CSRC)kdtree.o
	

OBJ=	$(COBJ) 

EXEC=	./bin/addr

$(EXEC):$(OBJ)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJ) $(LFLAGS)

$(COBJ) : %.o : %.c
	$(CC) $(INC) $(CFLAGS) -c $< -o $@

clean:
	rm $(COBJ)
	rm $(EXEC)
