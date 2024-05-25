#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#define ARR_LEN 8000000

int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double tstart, tend;
    
    int lg2, pow2;
    double chunk = ceil((double)ARR_LEN / size);
    int start = rank * chunk;
    int end = (rank+1) * chunk - 1;
    long* partial_sums = (long *)malloc(ARR_LEN * sizeof(long)); // Массив для хранения частных сумм
    long* sums = (long *)malloc(chunk * sizeof(long));
    

    lg2 = round(log2(ARR_LEN)) + 1;

    // Заполнение массива (просто для примера)
    for (int i = 0; i < ARR_LEN; ++i) 
        partial_sums[i] = i+1;

    if (rank == 0)
    {
        printf("input array\n");
        if (ARR_LEN > 20)
        {
            for (int i = 0; i < 20; ++i)
                printf("%ld ", partial_sums[i]);
            printf("...\n\n");
        }
        else
        {
            for (int i = 0; i < ARR_LEN; ++i)
                printf("%ld ", partial_sums[i]);
            printf("\n\n"); 
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    tstart = MPI_Wtime();

    for (int i = 0; i < lg2; i++)
    {
        pow2 = pow(2, i);
        
        for (int j = start; j <= end; j++)
        {
            if (j >= pow2)
                sums[j-start] = partial_sums[j] + partial_sums[j - pow2];
            else
                sums[j-start] = partial_sums[j];
        }
        MPI_Allgather(sums, chunk, MPI_LONG, partial_sums, chunk, MPI_LONG, MPI_COMM_WORLD);
    }

    tend = MPI_Wtime();

    if (rank == 0)
    {
        printf("partial_sums res\n");
        if (ARR_LEN > 20)
        {
            for (int i = 0; i < 20; ++i)
                printf("%ld ", partial_sums[i]);
            printf("...\n\n");
        }
        else
        {
            for (int i = 0; i < ARR_LEN; ++i)
                printf("%ld ", partial_sums[i]);
            printf("\n\n"); 
        }

        double time_taken_parallel = tend - tstart;
        printf("For array size %d\n", ARR_LEN);
        printf("Time spend parallel = %f\n", time_taken_parallel);
    }
    MPI_Finalize();
    return 0;
}