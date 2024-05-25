#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

struct pair {
    int first;
    int second;
};

struct pPair {
    pair first;
    double second;
};

struct cell {
    // Row and Column index of its parent
    // Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
    int parent_i, parent_j;
    // f = g + h
    double f, g, h;
};

int **allocate_matrix(int n, int m) {
    int **matrix = (int **)malloc(n * sizeof(int *));
    for (int i = 0; i < n; i++) {
        matrix[i] = (int *)malloc(m * sizeof(int));
    }
    return matrix;
}

void free_matrix(int **matrix, int n) {
    for (int i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

int readmatrix(size_t rows, size_t cols, int** matrix , const char* filename)
{

    FILE *pf;
    pf = fopen (filename, "r");
    if (pf == NULL)
        return 0;

    for(size_t i = 0; i < rows; ++i)
    {
        for(size_t j = 0; j < cols; ++j)
            fscanf(pf, "%d", matrix[i] + j);
    }


    fclose (pf); 
    return 1; 
}

int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int** path = allocate_matrix(5, 5);
    readmatrix(5, 5, path, "in_path.txt");

    MPI_Finalize();
    return 0;
}