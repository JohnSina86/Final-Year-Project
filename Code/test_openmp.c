#include <stdio.h>
#include <omp.h>

int main() {
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int num_threads = omp_get_num_threads();

        #pragma omp critical
        printf("Thread %d of %d\n", thread_id, num_threads);
    }

    printf("\nOpenMP is working!\n");
    return 0;
}