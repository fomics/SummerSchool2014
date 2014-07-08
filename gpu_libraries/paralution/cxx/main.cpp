// ******************************************
// implicit time stepping implementation of 2D diffusion problem
// Ben Cumming, CSCS
// *****************************************

// A small benchmark app that solves the 2D fisher equation using second-order
// finite differences.

// Syntax: ./main nx ny nt t

#include <algorithm>
#include <fstream>
#include <iostream>

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <omp.h>
#include <stdio.h>

#include <paralution.hpp>

#include "data.h"
#include "linalg.h"
#include "operators.h"
#include "matrix.h"
#include "stats.h"

using namespace data;
using namespace linalg;
using namespace operators;
using namespace matrix;
using namespace stats;
using namespace paralution;

// ==============================================================================

// read command line arguments
static void readcmdline(Discretization& options, int argc, char* argv[])
{
    if (argc<5 || argc>6)
    {
        printf("Usage: main nx ny nt t verbose\n");
        printf("  nx        number of gridpoints in x-direction\n");
        printf("  ny        number of gridpoints in y-direction\n");
        printf("  nt        number of timesteps\n");
        printf("  t         total time\n");
        printf("  verbose   (optional) verbose output\n");
        exit(1);
    }

    // read nx
    options.nx = atoi(argv[1]);
    if (options.nx < 1)
    {
        fprintf(stderr, "nx must be positive integer\n");
        exit(-1);
    }

    // read ny
    options.ny = atoi(argv[2]);
    if (options.ny < 1)
    {
        fprintf(stderr, "ny must be positive integer\n");
        exit(-1);
    }

    // read nt
    options.nt = atoi(argv[3]);
    if (options.nt < 1)
    {
        fprintf(stderr, "nt must be positive integer\n");
        exit(-1);
    }

    // read total time
    double t = atof(argv[4]);
    if (t < 0)
    {
        fprintf(stderr, "t must be positive real value\n");
        exit(-1);
    }

    // set verbosity if requested
    if (argc==6)
        verbose_output = true;
    else
        verbose_output = false;

    // store the parameters
    options.N = options.nx * options.ny;

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
    init_paralution();
    info_paralution();

    // read command line arguments
    readcmdline(options, argc, argv);
    int nx = options.nx;
    int ny = options.ny;
    int N  = options.N;
    int nt = options.nt;

    int unknowns = nx*ny;

    // Additions to form COO version of matrix
    int nnz = 5*nx*ny - 2*nx - 2*ny;
    LocalMatrix<double> mat ;
    LocalVector<double> rhs ;
    LocalVector<double> delta ;
    LocalVector<double> xold ;
    LocalVector<double> xnew ;

    int* Iindices       = new int[nnz];
    int* Jindices       = new int[nnz];
    double* Values      = new double[nnz];

    std::cout << "========================================================================" << std::endl;
    std::cout << "                      Welcome to mini-stencil!" << std::endl;
    std::cout << "version :: C++ serial" << std::endl;
    std::cout << "mesh    :: " << options.nx << " * " << options.ny << " dx = " << options.dx << std::endl;
    std::cout << "time    :: " << nt << " time steps from 0 .. " << options.nt*options.dt << std::endl;;
    std::cout << "========================================================================" << std::endl;

    // allocate global fields
    // allocate global fields
    x_new.init(nx,ny);
    x_old.init(nx,ny);
    bndN.init(nx,1);
    bndS.init(nx,1);
    bndE.init(ny,1);
    bndW.init(ny,1);

    Field b(nx,ny);
    Field deltax(nx,ny);

    // set dirichlet boundary conditions to 0 all around
    ss_fill(bndN, 0., nx);
    ss_fill(bndS, 0., nx);
    ss_fill(bndE, 0., ny);
    ss_fill(bndW, 0., ny);

    // set the initial condition
    // a circle of concentration 0.1 centred at (xdim/4, ydim/4) with radius
    // no larger than 1/8 of both xdim and ydim
    ss_fill(x_new, 0., unknowns);
    double xc = 1.0 / 4.0;
    double yc = (ny - 1) * options.dx / 4;
    double radius = fmin(xc, yc) / 2.0;
    for (int j = 0; j < ny; j++)
    {
        double y = (j - 1) * options.dx;
        for (int i = 0; i < nx; i++)
        {
            double x = (i - 1) * options.dx;
            if ((x - xc) * (x - xc) + (y - yc) * (y - yc) < radius * radius)
                x_new[i+nx*j] = 0.1;
        }
    }

    define_coo(Iindices, Jindices, Values);
    mat.Assemble( Iindices, Jindices, Values, nnz, "A" );
    xold.Assemble( Iindices, x_old.data(), unknowns, "rhs"); 

    mat.MoveToAccelerator();
    xold.MoveToAccelerator();
    rhs.MoveToAccelerator();
    delta.MoveToAccelerator();

    rhs.Allocate("rhs", unknowns);  // will be assembled subsequently
    delta.Allocate("rhs", unknowns);

    // Linear Solver
    CG<LocalMatrix<double>, LocalVector<double>, double > ls;

    ls.SetOperator(mat);
    ls.Build();
    
    iters_cg = 0;
    iters_newton = 0;

    // start timer
    double timespent = -omp_get_wtime();

    // main timeloop
    double tolerance = 1.e-6;
    for (int timestep = 1; timestep <= nt; timestep++)
    {
        // set x_new and x_old to be the solution
        ss_copy(x_old, x_new, N);
        xnew.CloneFrom( xold ); // PARALUTION

        double residual;
        bool converged = false;
        int it;
        for (it=0; it<50; it++)
        {
          // PARALUTION
	  delta.Zeros();
	  define_update_coo(x_new, Values);
          mat.AssembleUpdate(Values);   // Update the modified values 
          rhs.Assemble( Iindices, b.data(), unknowns, "rhs"); 
          ls.ResetOperator(mat);
          mat.info();

	  // compute residual : requires both x_new and x_old
	  diffusion(x_new, b);
	  // residual = ss_norm2(b, N);
          
	  // check for convergence
	  if (residual < tolerance && it > 0)
            {
                converged = true;
                break;
            }

	  // solve linear system to get -deltax
	  bool cg_converged = false;
	  ss_cg(deltax, b, 200, tolerance, cg_converged);
	  ls.Solve(rhs, &delta);

	  // check that the CG solver converged
	  if (!cg_converged) break;

	  // update solution
	  ss_axpy(x_new, -1.0, deltax, N);
          residual = delta.Norm();
          xnew.AddScale(delta, -1.0); // PARALUTION: add in delta

        }
        iters_newton += it+1;

        // output some statistics
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

    ////////////////////////////////////////////////////////////////////
    // write final solution to BOV file for visualization
    ////////////////////////////////////////////////////////////////////

    // binary data
    {
        FILE* output = fopen("output.bin", "w");
        fwrite(x_new.data(), sizeof(double), nx * ny, output);
        fclose(output);
    }

    std::ofstream fid("output.bov");
    fid << "TIME: 0.0" << std::endl;
    fid << "DATA_FILE: output.bin" << std::endl;
    fid << "DATA_SIZE: " << options.nx << ", " << options.ny << ", 1" << std::endl;;
    fid << "DATA_FORMAT: DOUBLE" << std::endl;
    fid << "VARIABLE: phi" << std::endl;
    fid << "DATA_ENDIAN: LITTLE" << std::endl;
    fid << "CENTERING: nodal" << std::endl;
    fid << "BRICK_SIZE: 1.0 " << (options.ny-1)*options.dx << " 1.0" << std::endl;

    // print table sumarizing results
    std::cout << "--------------------------------------------------------------------------------"
              << std::endl;
    std::cout << "simulation took " << timespent << " seconds" << std::endl;
    std::cout << int(iters_cg) << " conjugate gradient iterations, at rate of "
              << float(iters_cg)/timespent << " iters/second" << std::endl;
    std::cout << iters_newton << " newton iterations" << std::endl;
    std::cout << "--------------------------------------------------------------------------------"
              << std::endl;

    std::cout << "Goodbye!" << std::endl;

    stop_paralution();

    delete( Iindices );
    delete( Jindices );
    delete( Values   );

    return 0;
}

