FTN = gfortran
FFLAGS = -O3 -fopenmp
ifeq ($(PE_ENV),GNU)
	FTN=ftn
	FFLAGS=-O3 -fopenmp
endif
ifeq ($(PE_ENV),CRAY)
	FTN=ftn
	FFLAGS=-O3
endif
ifeq ($(PE_ENV),PGI)
	FTN=ftn
	FFLAGS=-O3 -mp
endif
ifeq ($(PE_ENV),INTEL)
	FTN=ftn
	FFLAGS=-O3 -openmp
endif

SOURCES = stats.f90 data.f90 operators.f90 linalg.f90
OBJ     = stats.o   data.o   operators.o   linalg.o

.SUFFIXES: .f90

all: main

stats.o: stats.f90
	$(FTN) $(FFLAGS) $(IFLAGS) -c stats.f90 -o stats.o

data.o: data.f90
	$(FTN) $(FFLAGS) $(IFLAGS) -c data.f90 -o data.o

operators.o: operators.f90
	$(FTN) $(FFLAGS) $(IFLAGS) -c operators.f90 -o operators.o

linalg.o: linalg.f90
	$(FTN) $(FFLAGS) $(IFLAGS) -c linalg.f90 -o linalg.o

main: $(OBJ) main.f90
	$(FTN) $(FFLAGS) *.o main.f90  -o $@

clean:
	rm -f main
	rm -f *.o
	rm -f *.i
	rm -f *.mod
	rm -f output.*

