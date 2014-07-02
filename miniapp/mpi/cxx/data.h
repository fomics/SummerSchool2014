#pragma once

#include <cassert>

namespace data
{

// define some helper types that can be used to pass simulation
// data around without haveing to pass individual parameters
struct Discretization
{
    int nx;       // x dimension
    int ny;       // y dimension
    int nt;       // number of time steps
    double dt;    // time step size
    double dx;    // distance between grid points
    double alpha; // dx^2/(D*dt)
};

struct SubDomain
{
    // initialize a subdomain
    void init(int rank, int size, Discretization& options);

    // print subdomain information
    void print();

    // i and j dimensions of the global decomposition
    int ndomx;
    int ndomy;

    // the i and j index of this sub-domain
    int domx;
    int domy;

    // the i and j bounding box of this sub-domain
    int startx;
    int starty;
    int endx;
    int endy;

    // the rank of neighbouring domains
    int neighbour_north;
    int neighbour_east;
    int neighbour_south;
    int neighbour_west;

    // mpi info
    int size;
    int rank;

    // x and y dimension in grid points of the sub-domain
    int nx;
    int ny;

    // total number of grid points
    int N;
};

// thin wrapper around a pointer that can be accessed as either a 2D or 1D array
// Field has dimension xdim * ydim in 2D, or length=xdim*ydim in 1D
class Field {
    public:
    // default constructor
    Field()
    :   ptr_(0), xdim_(0), ydim_(0)
    {};

    // constructor
    Field(int xdim, int ydim)
    :   xdim_(xdim), ydim_(ydim)
    {
        #ifdef DEBUG
        assert(xdim>0 && ydim>0);
        #endif
        ptr_ = new double[xdim*ydim];
        // do first touch
        fill(0.);
    };

    // destructor
    ~Field() {
        free();
    }

    void init(int xdim, int ydim) {
        #ifdef DEBUG
        assert(xdim>0 && ydim>0);
        #endif
        free();
        ptr_ = new double[xdim*ydim];
        xdim_ = xdim;
        ydim_ = ydim;

        // do first touch
        fill(0.);
    }

    double*       data()       { return ptr_; }
    const double* data() const { return ptr_; }

    // access via (i,j) pair
    inline double&       operator() (int i, int j)        {
        #ifdef DEBUG
        assert(i>=0 && i<xdim_ && j>=0 && j<ydim_);
        #endif
        return ptr_[i+j*xdim_];
    }
    inline double const& operator() (int i, int j) const  {
        #ifdef DEBUG
        assert(i>=0 && i<xdim_ && j>=0 && j<ydim_);
        #endif
        return ptr_[i+j*xdim_];
    }

    // access as a 1D field
    inline double      & operator[] (int i) {
        #ifdef DEBUG
        assert(i>=0 && i<xdim_*ydim_);
        #endif
        return ptr_[i];
    }
    inline double const& operator[] (int i) const {
        #ifdef DEBUG
        assert(i>=0 && i<xdim_*ydim_);
        #endif
        return ptr_[i];
    }

    int xdim()   const { return xdim_; }
    int ydim()   const { return ydim_; }
    int length() const { return xdim_*ydim_; }

    private:

    // set to a constant value
    void fill(double val) {
        #pragma omp parallel for
        for(int i=0; i<xdim_*ydim_; ++i)
            ptr_[i] = val;
    }

    void free() {
        if(ptr_) delete[] ptr_;
        ptr_ = 0;
    }

    double* ptr_;
    int xdim_;
    int ydim_;
};

// fields that hold the solution
extern Field x_new; // 2d
extern Field x_old; // 2d

// fields that hold the boundary values
extern Field bndN; // 1d
extern Field bndE; // 1d
extern Field bndS; // 1d
extern Field bndW; // 1d

// buffers used in boundary exchange
extern Field buffN;
extern Field buffE;
extern Field buffS;
extern Field buffW;

extern Discretization options;
extern SubDomain      domain;

} // namespace data

