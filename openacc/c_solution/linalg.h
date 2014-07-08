// linear algebra subroutines
// Ben Cumming @ CSCS

#ifndef LINALG_H
#define LINALG_H

#include "data.h"

    extern int cg_initialized;
    extern double *r, *Ap, *p, *Fx, *Fxold, *v, *xold; // 1d

    ////////////////////////////////////////////////////////////////////////////////
    //  blas level 1 reductions
    ////////////////////////////////////////////////////////////////////////////////

    // computes the inner product of x and y
    // x and y are vectors on length N
    double ss_dot(const double* const __restrict__ x, const double* const __restrict__ y, const int N);

    // computes the 2-norm of x
    // x is a vector on length N
    double ss_norm2(const double* const __restrict__ x, const int N);

    // sets entries in a vector to value
    // x is a vector on length N
    // value is th
    void ss_fill(double* __restrict__ x, const double value, const int N);

    ////////////////////////////////////////////////////////////////////////////////
    //  blas level 1 vector-vector operations
    ////////////////////////////////////////////////////////////////////////////////

    // computes y := alpha*x + y
    // x and y are vectors on length N
    // alpha is a scalar
    void ss_axpy(double* __restrict__ y, const double alpha, const double* const __restrict__ x, const int N);

    // computes y = x + alpha*(l-r)
    // y, x, l and r are vectors of length N
    // alpha is a scalar
    void ss_add_scaled_diff(double* __restrict__ y, const double* const __restrict__ x, const double alpha,
        const double* const __restrict__ l, const double* const __restrict__ r, const int N);

    // computes y = alpha*(l-r)
    // y, l and r are vectors of length N
    // alpha is a scalar
    void ss_scaled_diff(double* __restrict__ y, const double alpha,
        const double* const __restrict__ l, const double* const __restrict__ r, const int N);

    // computes y := alpha*x
    // alpha is scalar
    // y and x are vectors on length n
    void ss_scale(double* __restrict__ y, const double alpha, const double* const __restrict__ x, const int N);

    // computes linear combination of two vectors y := alpha*x + beta*z
    // alpha and beta are scalar
    // y, x and z are vectors on length n
    void ss_lcomb(double* __restrict__  y, const double alpha, const double* const __restrict__ x, const double beta,
        const double* const __restrict__ z, const int N);

    // copy one vector into another y := x
    // x and y are vectors of length N
    void ss_copy(double* __restrict__ y, const double* const __restrict__ x, const int N);

    // conjugate gradient solver
    // solve the linear system A*x = b for x
    // the matrix A is implicit in the objective function for the diffusion equation
    // the value in x constitute the "first guess" at the solution
    // x(N)
    // ON ENTRY contains the initial guess for the solution
    // ON EXIT  contains the solution
    void ss_cg(double* __restrict__ x, const double* const __restrict__ b, const int maxiters, const double tol, int* success);

#endif // LINALG_H

