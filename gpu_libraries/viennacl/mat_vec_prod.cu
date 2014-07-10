/* =========================================================================
   Copyright (c) 2010-2013, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at
               
   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


/*
*   Benchmark:  Sparse matrix operations, i.e. matrix-vector products; adapted for GPU-enabled Libraries Tutorial by W. Sawyer
*   
*/

//#define VIENNACL_BUILD_INFO
#ifndef NDEBUG
 #define NDEBUG
#endif

#define VIENNACL_WITH_UBLAS 1

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/numeric/ublas/lu.hpp>


// EXERCISE:  For your convenience, we have added all the include files you could possibly need for this exercise

#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/io/matrix_market.hpp"
#include "viennacl/linalg/ilu.hpp"


#include <iostream>
#include <vector>
#include "benchmark-utils.hpp"
#include "io.hpp"


#define BENCHMARK_RUNS          10

// EXERCISE: please note that this routine is templated with a type for the scalar, e.g. single or double precision

template<typename ScalarType>
int run_benchmark()
{   
   Timer timer;
   double exec_time;
   
   //ScalarType std_result = 0;
   
  boost::numeric::ublas::vector<ScalarType> ublas_vec1;
  boost::numeric::ublas::vector<ScalarType> ublas_vec2;

  if (!readVectorFromFile<ScalarType>("/project/csstaff/inputs/Matrices/MatrixMarket/AutumnSchool2013/result65025.txt", ublas_vec1))
  {
    std::cout << "Error reading RHS file" << std::endl;
    return 0;
  }
  std::cout << "done reading rhs" << std::endl;
  ublas_vec2 = ublas_vec1;

  // Define a ViennaCL compressed matrix (see User Manual) named "vcl_compressed_matrix"
  
  boost::numeric::ublas::compressed_matrix<ScalarType> ublas_matrix;
  // EXERCISE: We use ViennaCL utilities to read in the matrix into the ublas_matrix.
  if (!viennacl::io::read_matrix_market_file(ublas_matrix, "/project/csstaff/inputs/Matrices/MatrixMarket/AutumnSchool2013/mat65k.mtx"))
  {
    std::cout << "Error reading Matrix file" << std::endl;
    return 0;
  }
  std::cout << "done reading matrix" << std::endl;

  // EXERCISE:  here you need to define ViennaCL vectors vcl_vec1 and vcl_vec2, with size ublas_vec1.size()
  
  // EXERCISE:  copy cpu values to gpu

  //     a)  copy ublas_matrix to vcl_compressed_matrix
  //     b)  copy ublas_vec1   to vcl_vec1
  //     c)  copy ublas_vec2   to vcl_vec2

  ///////////// Matrix operations /////////////////
  
  std::cout << "------- Matrix-Vector product on CPU ----------" << std::endl;
  timer.start();
  for (int runs=0; runs<BENCHMARK_RUNS; ++runs)
  {
    //ublas_vec1 = boost::numeric::ublas::prod(ublas_matrix, ublas_vec2);
    boost::numeric::ublas::axpy_prod(ublas_matrix, ublas_vec2, ublas_vec1, true);
  }
  exec_time = timer.get();
  std::cout << "CPU time: " << exec_time << std::endl;
  std::cout << "CPU "; printOps(2.0 * static_cast<double>(ublas_matrix.nnz()), static_cast<double>(exec_time) / static_cast<double>(BENCHMARK_RUNS));
  std::cout << ublas_vec1[0] << std::endl;
  
  
  std::cout << "------- Matrix-Vector product with ViennaCL compressed_matrix ----------" << std::endl;
  
  // EXERCISE: here it is necessary to perform 1 mat-vec product to "warm up" the GPU
  //           vcl_vec1 = vcl_compressed_matrix * vcl_vec2
  
  viennacl::backend::finish();   // This completes the viennacl execution

  timer.start();
  for (int runs=0; runs<BENCHMARK_RUNS; ++runs)
  {
    // EXERCISE: perform mat-vec product
    //           vcl_vec1 = vcl_compressed_matrix * vcl_vec2
    std::cout << "===>  Replace this output line with the matrix-vector product <===" << std::endl;
  }
  viennacl::backend::finish();
  exec_time = timer.get();
  std::cout << "GPU time: " << exec_time << std::endl;
  std::cout << "GPU "; printOps(2.0 * static_cast<double>(ublas_matrix.nnz()), static_cast<double>(exec_time) / static_cast<double>(BENCHMARK_RUNS));

  return EXIT_SUCCESS;
}


int main()
{
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "## Benchmark :: Sparse" << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << std::endl;
  std::cout << "   -------------------------------" << std::endl;
  std::cout << "   # benchmarking single-precision" << std::endl;
  std::cout << "   -------------------------------" << std::endl;
  run_benchmark<float>();
  std::cout << std::endl;
  std::cout << "   -------------------------------" << std::endl;
  std::cout << "   # benchmarking double-precision" << std::endl;
  std::cout << "   -------------------------------" << std::endl;
  run_benchmark<double>();
  return 0;
}
