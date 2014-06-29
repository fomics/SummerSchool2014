#include <iostream>

#include <omp.h>

int main(void){
    int tid;

    std::cout << "=== serial section ===" << std::endl;

    std::cout << "hello world from thread " << omp_get_thread_num()
              << " of " << omp_get_num_threads()
              << std::endl;

    std::cout << "=== parallel section ===" << std::endl;
    #pragma omp parallel
    {
        // get the number of this thread
        tid = omp_get_thread_num();

        // write a personalized message from this thread
        std::cout << "hello world from thread " << tid
                  << " of " << omp_get_num_threads()
                  << std::endl;
    }

    return 0;
}
