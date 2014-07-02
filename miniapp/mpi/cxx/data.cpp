#include <iostream>

#include <cmath>

#include <mpi.h>
#include <omp.h>

#include "data.h"

namespace data{

// fields that hold the solution
Field x_new;
Field x_old;

// fields that hold the boundary points
Field bndN;
Field bndE;
Field bndS;
Field bndW;

// buffers used during boundary halo communication
Field buffN;
Field buffE;
Field buffS;
Field buffW;

Discretization options;
SubDomain      domain;

void SubDomain::init(int mpi_rank, int mpi_size, Discretization& discretization)
{
    // determine the number of subdomains in the x and y dimensions
    ndomx = sqrt(double(mpi_size));
    while( mpi_size%ndomx )
        ndomx--;
    ndomy = mpi_size / ndomx;

    // compute this sub-domain index
    // work backwards from: mpi_rank = (domx-1) + (domy-1)*ndomx
    domx = mpi_rank % ndomx + 1;
    domy = (mpi_rank-domx+1) / ndomx + 1;

    nx = discretization.nx / ndomx;
    ny = discretization.ny / ndomy;
    startx = (domx-1)*nx+1;
    starty = (domy-1)*ny+1;

    // adjust for grid dimensions that do not divided evenly between the
    // sub-domains
    if( domx == ndomx )
        nx = discretization.nx - startx + 1;
    if( domy == ndomy )
        ny = discretization.ny - starty + 1;

    endx = startx + nx -1;
    endy = starty + ny -1;

    // get total number of grid points in this sub-domain
    N = nx*ny;

    ////////////////////////////////////////////////////////////////////
    // TODOSS
    // determine the ranks of your 4 neighbours
    ////////////////////////////////////////////////////////////////////
    neighbour_east  = -2;
    neighbour_west  = -2;
    neighbour_north = -2;
    neighbour_south = -2;
    
    // TODOSS
    // set rank of neighbours that lie outside the domain to be -1
    // neighbour_west, neighbour_east
    // hint: you will need domx
    if (domy == 1) {
        neighbour_south = -1;
    }
    if (domy == ndomy) {
        neighbour_north = -1;
    }

    rank = mpi_rank;
    size = mpi_size;
}

// print domain decomposition information to stdout
void SubDomain::print() {
    if(rank==0) {
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::cout << "there are " << size << " subdomains on a "
                  << ndomx << " * " << ndomy << " grid " << std::endl;
        std::cout << "global grid has dimensions "
                  <<  options.nx << " * " << options.ny << std::endl;
        std::cout << "-----------------------------------------------------------------------" << std::endl;
    }
    for(int i=0; i<size; i++) {
        if(rank == i) {
            std::cout << "rank " << rank
                      << " (" << domx << "," << domy << ")"
                      << " : neighbours N:S:E:W (" << neighbour_north << "  " << neighbour_south
                      << "  " << neighbour_east << "  " << neighbour_west
                      << ") : local dims " << nx << " * " << ny
                      << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    double time = omp_get_wtime();
    // add artificial pause so that output doesn't pollute later output
    while(omp_get_wtime()-time < 1e-1);
}

} // namespace data
