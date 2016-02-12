###################################################################
#
# freeman.justin@gmail.com 
#
##################################################################

CC=	gcc

CSRC=	./src/

CFLAGS=	-O3 -g -Wall `xml2-config --cflags`
#CFLAGS=	-O3 -g -fPIC -Wall

NNLIB=	./lib/nn

INC=	-I./include \
	-I$(NNLIB)


#LFLAGS= -lnetcdf -lxml2 
LFLAGS= -lxml2 -ludunits2 -lnetcdf $(NNLIB)/libnn.a

COBJ=	$(CSRC)main.o \
	$(CSRC)jutil.o \
	$(CSRC)fail.o \
	$(CSRC)spheriq_dist.o \
	$(CSRC)malloc_arrays.o \
	$(CSRC)netcdfIO.o \
	$(CSRC)interp.o \
	$(CSRC)calcalcs.o \
	$(CSRC)utCalendar2_cal.o \
	$(CSRC)nnutil.o \
	$(CSRC)process_roms.o \
	$(CSRC)process_tides.o \
	$(CSRC)process_waves.o \
	$(CSRC)add.o
	

OBJ=	$(COBJ) 

EXEC=	./bin/addr

$(EXEC):$(OBJ)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJ) $(LFLAGS)

$(COBJ) : %.o : %.c
	$(CC) $(INC) $(CFLAGS) -c $< -o $@

clean:
	rm $(COBJ)
	rm $(EXEC)
