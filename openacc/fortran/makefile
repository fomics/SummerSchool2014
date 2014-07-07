FTN = ftn

FFLAGS =
ifeq ($(PE_ENV),PGI)
	# -Mifno outputs extra diagnostic information
	# -hacc generates accelerated kernels from openacc directives
	FFLAGS  = -acc=verystrict
	FFLAGS += -ta=nvidia,nofma,cc35,cuda5.5 # target architecture
	#FFLAGS  = -Minfo=accel
	FFLAGS += -Mpreprocess # use c-style preprocessing
	FFLAGS += -g # debug flags
endif
ifeq ($(PE_ENV),CRAY)
	# -rmd outputs extra diagnostic files
	# -hacc generates accelerated kernels from openacc directives
	FFLAGS  = -rmd -hacc -O3
	FFLAGS += -e Z # use c-style preprocessing
endif

SOURCES = stats.f90 data.f90 operators.f90 linalg.f90
OBJ     = stats.o   data.o   operators.o     linalg.o     io.o
OBJ_MPI = stats.o   data.o   operators_mpi.o linalg.o     io.o

.SUFFIXES: .f90

all: main main_mpi

stats.o: stats.f90
	$(FTN) $(FFLAGS) $(IFLAGS) -c stats.f90 -o stats.o

data.o: data.f90
	$(FTN) $(FFLAGS) $(IFLAGS) -c data.f90 -o data.o

io.o: io.f90
	$(FTN) $(FFLAGS) $(IFLAGS) -c io.f90 -o io.o

operators.o: operators.f90
	$(FTN) $(FFLAGS) $(IFLAGS) -c operators.f90 -o operators.o

operators_mpi.o: operators_mpi.f90
	$(FTN) $(FFLAGS) $(IFLAGS) -c operators_mpi.f90 -o operators_mpi.o
#	$(FTN) $(FFLAGS) $(IFLAGS) -c operators_mpi.f90 -o operators_mpi.o

linalg.o: linalg.f90
	$(FTN) $(FFLAGS) $(IFLAGS) -c linalg.f90 -o linalg.o

main: $(OBJ) main.f90
	$(FTN) $(FFLAGS) $(OBJ) main.f90  -o main

main_mpi: $(OBJ_MPI) main.f90
	$(FTN) $(FFLAGS) $(OBJ_MPI) main.f90  -o main_mpi

clean:
	rm -f main
	rm -f *.o
	rm -f *.i
	rm -f *.mod
	rm -f *.lst
	rm -f *.opt
	rm -f *.cg
	rm -f *.ptx
	rm -f *.cub
	rm -f output.*

