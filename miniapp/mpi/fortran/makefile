FTN = ftn

# set compiler-specific flags
ifeq ($(PE_ENV),GNU)
# the gnu compiler requires that we specifically enable OpenMP
	FFLAGS=-O2 -dynamic -fopenmp
endif
ifeq ($(PE_ENV),CRAY)
#	FFLAGS= -O2 -dynamic
endif


SOURCES = stats.f90 data.f90 operators.f90 linalg.f90
OBJ     = stats.o   data.o   operators.o   linalg.o     io.o

.SUFFIXES: .f90

all: main

stats.o: stats.f90
	$(FTN) $(FFLAGS) $(IFLAGS) -c stats.f90 -o stats.o

data.o: data.f90
	$(FTN) $(FFLAGS) $(IFLAGS) -c data.f90 -o data.o

io.o: io.f90
	$(FTN) $(FFLAGS) $(IFLAGS) -c io.f90 -o io.o

operators.o: operators.f90
	$(FTN) $(FFLAGS) $(IFLAGS) -c operators.f90 -o operators.o

linalg.o: linalg.f90
	$(FTN) $(FFLAGS) $(IFLAGS) -c linalg.f90 -o linalg.o

main: $(OBJ) main.f90
	$(FTN) $(FFLAGS) *.o main.f90  -o main
#	$(FTN) $(FFLAGS) *.o main.f90  -o $@ /apps/santis/ddt/4.2.1-35513/lib/64/libdmallocth.a

clean:
	rm -f main
	rm -f *.o
	rm -f *.i
	rm -f *.mod
	rm -f output.*

