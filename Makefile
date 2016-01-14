###################################################################
#
# freeman.justin@gmail.com 
#
##################################################################

CC=	gcc

CSRC=	./src/

CFLAGS=	-O3 -g -Wall `xml2-config --cflags`
#CFLAGS=	-O3 -g -fPIC -Wall

INC=	-I./include 

#LFLAGS= -lnetcdf -lxml2 
LFLAGS= -lxml2 -ludunits2 -lnetcdf

COBJ=	$(CSRC)main.o \
	$(CSRC)jutil.o \
	$(CSRC)fail.o \
	$(CSRC)spheriq_dist.o \
	$(CSRC)malloc_arrays.o \
	$(CSRC)netcdfIO.o \
	$(CSRC)xmlIO.o \
	$(CSRC)interp.o \
	$(CSRC)calcalcs.o \
	$(CSRC)utCalendar2_cal.o
	

OBJ=	$(COBJ) 

EXEC=	./bin/addr

$(EXEC):$(OBJ)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJ) $(LFLAGS)

$(COBJ) : %.o : %.c
	$(CC) $(INC) $(CFLAGS) -c $< -o $@

clean:
	rm $(COBJ)
	rm $(EXEC)
