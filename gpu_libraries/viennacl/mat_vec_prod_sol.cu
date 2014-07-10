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
*   Benchmark:  Sparse matrix operations, i.e. matrix-vector products (sparse.cpp and sparse.cu are identical, the latter being required for compilation using CUDA nvcc)
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


template<typename ScalarType>
int run_benchmark()
{   
   Timer timer;
   double exec_time;
   
   //ScalarType std_result = 0;
   
  boost::numeric::ublas::vector<ScalarType> ublas_vec1;
  boost::numeric::ublas::vector<ScalarType> ublas_vec2;

  if (!readVectorFromFile<ScalarType>("../testdata/result65025.txt", ublas_vec1))
  {
    std::cout << "Error reading RHS file" << std::endl;
    return 0;
  }
  std::cout << "done reading rhs" << std::endl;
  ublas_vec2 = ublas_vec1;
  
  viennacl::compressed_matrix<ScalarType, 1> vcl_compressed_matrix_1;
  
  boost::numeric::ublas::compressed_matrix<ScalarType> ublas_matrix;
  if (!viennacl::io::read_matrix_market_file(ublas_matrix, "../testdata/mat65k.mtx"))
  {
    std::cout << "Error reading Matrix file" << std::endl;
    return 0;
  }
  //unsigned int cg_mat_size = cg_mat.size(); 
  std::cout << "done reading matrix" << std::endl;
  
  viennacl::vector<ScalarType> vcl_vec1(ublas_vec1.size());
  viennacl::vector<ScalarType> vcl_vec2(ublas_vec1.size()); 
  
  //cpu to gpu:
  viennacl::copy(ublas_matrix, vcl_compressed_matrix_1);
  viennacl::copy(ublas_vec1, vcl_vec1);
  viennacl::copy(ublas_vec2, vcl_vec2);

  
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
  
  
  std::cout << "------- Matrix-Vector product with compressed_matrix ----------" << std::endl;
  
  
  vcl_vec1 = viennacl::linalg::prod(vcl_compressed_matrix_1, vcl_vec2); //startup calculation
  //std_result = 0.0;
  
  viennacl::backend::finish();
  timer.start();
  for (int runs=0; runs<BENCHMARK_RUNS; ++runs)
  {
    vcl_vec1 = viennacl::linalg::prod(vcl_compressed_matrix_1, vcl_vec2);
  }
  viennacl::backend::finish();
  exec_time = timer.get();
  std::cout << "GPU time align1: " << exec_time << std::endl;
  std::cout << "GPU align1 "; printOps(2.0 * static_cast<double>(ublas_matrix.nnz()), static_cast<double>(exec_time) / static_cast<double>(BENCHMARK_RUNS));
  std::cout << vcl_vec1[0] << std::endl;

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

