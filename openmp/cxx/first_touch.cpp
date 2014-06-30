#include <omp.h>
#include <cstdlib>
#include <iostream>

int main(void) {
    const size_t N = 100000000;
    double *x = (double*)malloc(sizeof(double)*N);
    double *y = (double*)malloc(sizeof(double)*N);
    double *z = (double*)malloc(sizeof(double)*N);

    for(int i=0; i<N; ++i) {
        x[i] = i;
        y[i] = i;
        z[i] = i;
    }

    double time = -omp_get_wtime();

    #pragma omp parallel for
    for(int i=0; i<N; ++i)
        z[i] += x[i]*y[i];

    time += omp_get_wtime();

    double sum = 0.;
    for(int i=0; i<N; ++i)
        sum += z[i];

    std::cout << sum << " took " << time << " seconds to compute" << std::endl;

}
