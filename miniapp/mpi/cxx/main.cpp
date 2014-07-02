// ******************************************
// implicit time stepping implementation of 2D diffusion problem
// Ben Cumming, CSCS
// *****************************************

// A small benchmark app that solves the 2D fisher equation using second-order
// finite differences.

// Syntax: ./main nx ny nt t
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <omp.h>
#include <mpi.h>

#include "data.h"
#include "linalg.h"
#include "operators.h"
#include "stats.h"

using namespace data;
using namespace linalg;
using namespace operators;
using namespace stats;

// ==============================================================================
void write_binary(std::string fname, Field &u, SubDomain &domain, Discretization &options)
{
    MPI_Offset disp = 0;
    MPI_File filehandle;
    MPI_Datatype filetype;

    int result =
        MPI_File_open(
            MPI_COMM_WORLD,
            fname.c_str(),
            MPI_MODE_CREATE | MPI_MODE_WRONLY,
            MPI_INFO_NULL,
            &filehandle
        );
    assert(result==MPI_SUCCESS);

    int ustart[]  = {domain.startx-1, domain.starty-1};
    int ucount[]  = {domain.nx, domain.ny};
    int dimuids[] = {options.nx, options.ny};

    result = MPI_Type_create_subarray(2, dimuids, ucount, ustart, MPI_ORDER_FORTRAN, MPI_DOUBLE, &filetype);
    assert(result==MPI_SUCCESS);

    result = MPI_Type_commit(&filetype);
    assert(result==MPI_SUCCESS);

    result = MPI_File_set_view(filehandle, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    assert(result==MPI_SUCCESS);

    result = MPI_File_write_all(filehandle, u.data(), domain.N, MPI_DOUBLE, MPI_STATUS_IGNORE);
    assert(result==MPI_SUCCESS);

    result = MPI_Type_free(&filetype);
    assert(result==MPI_SUCCESS);

    result = MPI_File_close(&filehandle);
    assert(result==MPI_SUCCESS);
}

// read command line arguments
static void readcmdline(Discretization& options, int argc, char* argv[])
{
    if (argc<5 || argc>6 ) {
        std::cerr << "Usage: main nx ny nt t\n";
        std::cerr << "  nx  number of gridpoints in x-direction\n";
        std::cerr << "  ny  number of gridpoints in y-direction\n";
        std::cerr << "  nt  number of timesteps\n";
        std::cerr << "  t   total time\n";
        std::cerr << "  v   [optional] turn on verbose output\n";
        exit(1);
    }

    // read nx
    options.nx = atoi(argv[1]);
    if (options.nx < 1) {
        std::cerr << "nx must be positive integer\n";
        exit(-1);
    }

    // read ny
    options.ny = atoi(argv[2]);
    if (options.ny < 1) {
        std::cerr << "ny must be positive integer\n";
        exit(-1);
    }

    // read nt
    options.nt = atoi(argv[3]);
    if (options.nt < 1) {
        std::cerr << "nt must be positive integer\n";
        exit(-1);
    }

    // read total time
    double t = atof(argv[4]);
    if (t < 0) {
        std::cerr << "t must be positive real value\n";
        exit(-1);
    }

    if( argc==6 )
        verbose_output = true;

    // compute timestep size
    options.dt = t / options.nt;

    // compute the distance between grid points
    // assume that x dimension has length 1.0
    options.dx = 1. / (options.nx - 1);

    // set alpha, assume diffusion coefficient D is 1
    options.alpha = (options.dx * options.dx) / (1. * options.dt);
}

// ==============================================================================

int main(int argc, char* argv[])
{
    // read command line arguments
    readcmdline(options, argc, argv);

    // initialize MPI
    int mpi_rank, mpi_size;
    if( MPI_Init(&argc, &argv) != MPI_SUCCESS ) {
        std::cerr << "unable to initialize MPI :: exitting" << std::endl;
        exit(-1);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // initialize subdomain
    domain.init(mpi_rank, mpi_size, options);
    if(verbose_output)
        domain.print();

    int nx = domain.nx;
    int ny = domain.ny;
    int N  = domain.N;
    int nt  = options.nt;

    if( domain.rank == 0 ) {
        std::cout << "========================================================================" << std::endl;
        std::cout << "                      Welcome to mini-stencil!" << std::endl;
        std::cout << "version :: C++ with MPI : " << domain.size << " MPI ranks" << std::endl;
        std::cout << "mesh    :: " << options.nx << " * " << options.ny << " dx = " << options.dx << std::endl;
        std::cout << "time    :: " << nt << " time steps from 0 .. " << options.nt*options.dt << std::endl;;
        std::cout << "========================================================================" << std::endl;
    }

    // allocate global fields
    x_new.init(nx,ny);
    x_old.init(nx,ny);
    bndN.init(nx,1);
    bndS.init(nx,1);
    bndE.init(ny,1);
    bndW.init(ny,1);
    buffN.init(nx,1);
    buffS.init(nx,1);
    buffE.init(ny,1);
    buffW.init(ny,1);

    // set dirichlet boundary conditions to 0 all around
    ss_fill(bndN, 0., nx);
    ss_fill(bndS, 0., nx);
    ss_fill(bndE, 0., ny);
    ss_fill(bndW, 0., ny);


    Field b(nx,ny);
    Field deltax(nx,ny);

    // set the initial condition
    // a circle of concentration 0.1 centred at (xdim/4, ydim/4) with radius
    // no larger than 1/8 of both xdim and ydim
    ss_fill(x_new, 0., nx*ny);
    double xc = 1.0 / 4.0;
    double yc = (options.ny - 1) * options.dx / 4;
    double radius = std::min(xc, yc) / 2.0;
    for (int j = domain.starty-1; j < domain.endy; j++)
    {
        double y = (j - 1) * options.dx;
        for (int i = domain.startx-1; i < domain.endx; i++)
        {
            double x = (i - 1) * options.dx;
            if ((x - xc) * (x - xc) + (y - yc) * (y - yc) < radius * radius)
                x_new(i-domain.startx+1, j-domain.starty+1) = 0.1;
        }
    }

    double time_in_bcs = 0.0;
    double time_in_diff = 0.0;
    flops_bc = 0;
    flops_diff = 0;
    flops_blas1 = 0;
    verbose_output = false;
    iters_cg = 0;
    iters_newton = 0;

    // start timer
    double timespent = -omp_get_wtime();

    // main timeloop
    double alpha = options.alpha;
    double tolerance = 1.e-6;
    for (int timestep = 1; timestep <= nt; timestep++)
    {
        // set x_new and x_old to be the solution
        ss_copy(x_old, x_new, N);

        double residual;
        bool converged = false;
        int it;
        for (it=0; it<50; it++)
        {
            // compute residual : requires both x_new and x_old
            diffusion(x_new, b);
            residual = ss_norm2(b, N);

            // check for convergence
            if (residual < tolerance)
            {
                converged = true;
                break;
            }

            // solve linear system to get -deltax
            bool cg_converged = false;
            ss_cg(deltax, b, 200, tolerance, cg_converged);

            // check that the CG solver converged
            if (!cg_converged) break;

            // update solution
            ss_axpy(x_new, -1.0, deltax, N);
        }
        iters_newton += it+1;

        #ifdef OUTPUT_EVERY_STEP
        std::stringstream str;
        str << "output" << timestep << ".bin"; // get filename for this time step
        write_binary(str.str(), x_old, domain, options); // write binary solution to file
        #endif

        // output some statistics
        //if (converged && verbose_output)
        if (converged && verbose_output) {
            std::cout << "step " << timestep
                      << " required " << it
                      << " iterations for residual " << residual
                      << std::endl;
        }
        if (!converged) {
            std::cerr << "step " << timestep
                      << " ERROR : nonlinear iterations failed to converge" << std::endl;;
            break;
        }
    }

    // get times
    timespent += omp_get_wtime();
    unsigned long long flops_total = flops_diff + flops_blas1;

    ////////////////////////////////////////////////////////////////////
    // write final solution to BOV file for visualization
    ////////////////////////////////////////////////////////////////////

    // binary data
    write_binary("output.bin", x_old, domain, options);

    // metadata
    if( domain.rank==0 ) {
        std::ofstream fid("output.bov");
        fid << "TIME: 0.0" << std::endl;
        fid << "DATA_FILE: output.bin" << std::endl;
        fid << "DATA_SIZE: " << options.nx << ", " << options.ny << ", 1" << std::endl;;
        fid << "DATA_FORMAT: DOUBLE" << std::endl;
        fid << "VARIABLE: phi" << std::endl;
        fid << "DATA_ENDIAN: LITTLE" << std::endl;
        fid << "CENTERING: nodal" << std::endl;
        fid << "BRICK_SIZE: 1.0 " << (options.ny-1)*options.dx << " 1.0" << std::endl;
    }

    // print table sumarizing results
    if(domain.rank == 0) {
        std::cout << "--------------------------------------------------------------------------------"
                  << std::endl;
        std::cout << "simulation took " << timespent << " seconds" << std::endl;
        std::cout << int(iters_cg) << " conjugate gradient iterations, at rate of "
                  << float(iters_cg)/timespent << " iters/second" << std::endl;
        std::cout << iters_newton << " newton iterations" << std::endl;
        std::cout << "--------------------------------------------------------------------------------"
                  << std::endl;
    }

    if(domain.rank==0)
        std::cout << "Goodbye!" << std::endl;

    MPI_Finalize();

    return 0;
}

