#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#define ROW 2048
#define COL 2048
#define NUM_SEARCH 4

typedef struct {
    int first, second;
}Pair;

Pair newPair(int first, int second) {
    Pair p;
    p.first = first;
    p.second = second;
    return p;
}

typedef struct {
    int row, col;
    double f;
}OpenCellData;

OpenCellData newOpenCellData(int row, int col, double f) {
    OpenCellData data;
    data.row = row;
    data.col = col;
    data.f = f;
    return data;
}

typedef struct{
    int parent_i, parent_j;
    double f, g, h;
}cell;

cell newCell(int parent_i, int parent_j, double f, double g, double h) {
    cell c;
    c.parent_i = parent_i;
    c.parent_j = parent_j;
    c.f = f;
    c.g = g;
    c.h = h;
    return c;
}

bool isValid(int row, int col)
{
    // Returns true if row number and column number is in range
    return (row >= 0) && (row < ROW) && (col >= 0) && (col < COL);
}
 
// A Utility Function to check whether the given cell is blocked or not
bool isUnBlocked(int** grid, int row, int col)
{
    // Returns true if the cell is not blocked else false
    if (grid[row][col] == 1)
        return (true);
    else
        return (false);
}
 
// A Utility Function to check whether destination cell has been reached or not
bool isDestination(int row, int col, int dest_row, int dest_col)
{
    if (row == dest_row && col == dest_col)
        return (true);
    else
        return (false);
}
 
// A Utility Function to calculate the 'h' heuristics.
double calculateHValue(int row, int col, int dest_row, int dest_col)
{
    // Return using the distance formula
    return ((double)sqrt(
        (row - dest_row) * (row - dest_row)
        + (col - dest_col) * (col - dest_col)));
}

void tracePath(cell** cellDetails, int dest_row, int dest_col, int** grid)
{
    int row = dest_row;
    int col = dest_col;
    Pair* Path = (Pair *)malloc(ROW * COL * sizeof(Pair));
    int size = 0;
    while (!(cellDetails[row][col].parent_i == row && cellDetails[row][col].parent_j == col)) 
    {
        size += 1;
        for (int i = size - 1; i > 0; i--)
        {
            Path[i] = Path[i - 1];
        }
        Path[0] = newPair(row, col);;
        
        int temp_row = cellDetails[row][col].parent_i;
        int temp_col = cellDetails[row][col].parent_j;
        row = temp_row;
        col = temp_col;
    }

    size += 1;
    for (int i = size - 1; i > 0; i--)
    {
        Path[i] = Path[i - 1];
    }
    Path[0] = newPair(row, col);

    for (int i = 0; i < size; i++)
    {
        grid[Path[i].first][Path[i].second] = 2;
    }
    return;
}

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

int writematrix(size_t rows, size_t cols, int** matrix , const char* filename)
{

    FILE *pf;
    pf = fopen (filename, "w+");
    if (pf == NULL)
        return 0;

    for(size_t i = 0; i < rows; ++i)
    {
        for(size_t j = 0; j < cols; ++j)
            fprintf(pf, "%d ", matrix[i][j]);
        fprintf(pf, "\n");
    }

    fclose (pf); 
    return 1; 
}

void aStar(int** grid, int** out_grid, int src_row, int src_col, int dest_row, int dest_col)
{
    // If the source is out of range
    if (isValid(src_row, src_col) == false) {
        printf("Source is invalid\n");
        return;
    }
 
    // If the destination is out of range
    if (isValid(dest_row, dest_col) == false) {
        printf("Destination is invalid\n");
        return;
    }
 
    // Either the source or the destination is blocked
    if (isUnBlocked(grid, src_row, src_col) == false
        || isUnBlocked(grid, dest_row, dest_col)
               == false) {
        printf("Source or the destination is blocked\n");
        return;
    }
 
    // If the destination cell is the same as source cell
    if (isDestination(src_row, src_col, dest_row, dest_col)
        == true) {
        printf("We are already at the destination\n");
        return;
    }
 
    // Create a closed list and initialise it to false which
    // means that no cell has been included yet This closed
    // list is implemented as a boolean 2D array

    bool closedList[ROW][COL];
    memset(closedList, false, sizeof(closedList));
 
    // Declare a 2D array of structure to hold the details
    // of that cell
    cell** cellDetails = (cell **)malloc(ROW * sizeof(cell *));
    for (int i = 0; i < ROW; ++i)
        cellDetails[i] = (cell *)malloc(COL * sizeof(cell));
        
    int i, j;
    for (i = 0; i < ROW; i++) {
        for (j = 0; j < COL; j++) {
            cellDetails[i][j] = newCell(-1, -1, FLT_MAX, FLT_MAX, FLT_MAX);
        }
    }

    // Initialising the parameters of the starting node
    i = src_row, j = src_col;
    cellDetails[i][j] = newCell(i, j, 0.0, 0.0, 0.0);

    // Put the starting cell on the open list and set its
    // 'f' as 0+
    OpenCellData* openList = (OpenCellData *)malloc(ROW * COL * sizeof(OpenCellData));
    openList[0] = newOpenCellData(i, j, 0.0);
    int openListSize = 1;
    // We set this boolean value as false as initially
    // the destination is not reached.
    bool foundDest = false;

    while (openListSize > 0) {
        OpenCellData p = openList[0];
 
        // Remove this vertex from the open list
        openListSize -= 1;
        for(i = 0; i < openListSize; i++)
            openList[i] = openList[i + 1];
        
        // Add this vertex to the closed list
        i = p.row;
        j = p.col;
        closedList[i][j] = true;
 
        // To store the 'g', 'h' and 'f' of the 8 successors
        double gNew, hNew, fNew;
 
        //----------- 1st Successor (North) ------------
        if (isValid(i - 1, j)) 
        {
            if (isDestination(i - 1, j, dest_row, dest_col) == true) 
            {
                cellDetails[i - 1][j].parent_i = i;
                cellDetails[i - 1][j].parent_j = j;
                printf("The destination cell is found\n");
                tracePath(cellDetails, dest_row, dest_col, out_grid);
                foundDest = true;
                return;
            }
            else if (closedList[i - 1][j] == false && isUnBlocked(grid, i - 1, j) == true)
            {
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue(i - 1, j, dest_row, dest_col);
                fNew = gNew + hNew;

                if (cellDetails[i - 1][j].f == FLT_MAX || cellDetails[i - 1][j].f > fNew) 
                {
                    openListSize += 1;
                    openList[openListSize - 1] = newOpenCellData(i - 1, j, fNew);
                    cellDetails[i - 1][j] = newCell(i, j, fNew, gNew, hNew);
                }
            }
        }

        //----------- 2nd Successor (South) ------------
        if (isValid(i + 1, j) == true) {
            if (isDestination(i + 1, j, dest_row, dest_col) == true)
            {
                cellDetails[i + 1][j].parent_i = i;
                cellDetails[i + 1][j].parent_j = j;
                printf("The destination cell is found\n");
                tracePath(cellDetails, dest_row, dest_col, out_grid);
                foundDest = true;
                return;
            }
            else if (closedList[i + 1][j] == false && isUnBlocked(grid, i + 1, j) == true)
            {
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue(i + 1, j, dest_row, dest_col);
                fNew = gNew + hNew;
                if (cellDetails[i + 1][j].f == FLT_MAX || cellDetails[i + 1][j].f > fNew)
                {
                    openListSize += 1;
                    openList[openListSize - 1] = newOpenCellData(i + 1, j, fNew);
                    cellDetails[i + 1][j] = newCell(i, j, fNew, gNew, hNew);
                }
            }
        }

        //----------- 3rd Successor (East) ------------
        if (isValid(i, j + 1) == true) 
        {
            if (isDestination(i, j + 1, dest_row, dest_col) == true) 
            {
                cellDetails[i][j + 1].parent_i = i;
                cellDetails[i][j + 1].parent_j = j;
                printf("The destination cell is found\n");
                tracePath(cellDetails, dest_row, dest_col, out_grid);
                foundDest = true;
                return;
            }
            else if (closedList[i][j + 1] == false && isUnBlocked(grid, i, j + 1) == true)
            {
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue(i, j + 1, dest_row, dest_col);
                fNew = gNew + hNew;
                if (cellDetails[i][j + 1].f == FLT_MAX || cellDetails[i][j + 1].f > fNew) 
                {
                    openListSize += 1;
                    openList[openListSize - 1] = newOpenCellData(i, j + 1, fNew);
                    cellDetails[i][j + 1] = newCell(i, j, fNew, gNew, hNew);
                }
            }
        }
        //----------- 4th Successor (West) ------------
        if (isValid(i, j - 1) == true) 
        {
            if (isDestination(i, j - 1, dest_row, dest_col) == true) {
                cellDetails[i][j - 1].parent_i = i;
                cellDetails[i][j - 1].parent_j = j;
                printf("The destination cell is found\n");
                tracePath(cellDetails, dest_row, dest_col, out_grid);
                foundDest = true;
                return;
            }
            else if (closedList[i][j - 1] == false && isUnBlocked(grid, i, j - 1) == true)
            {
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue(i, j - 1, dest_row, dest_col);
                fNew = gNew + hNew;
                if (cellDetails[i][j - 1].f == FLT_MAX || cellDetails[i][j - 1].f > fNew) 
                {
                    openListSize += 1;
                    openList[openListSize - 1] = newOpenCellData(i, j - 1, fNew);
                    cellDetails[i][j - 1] = newCell(i, j, fNew, gNew, hNew);
                }
            }
        }
        //----------- 5th Successor (North-East) ------------ 
        if (isValid(i - 1, j + 1) == true)
        {
            if (isDestination(i - 1, j + 1, dest_row, dest_col) == true) 
            {
                cellDetails[i - 1][j + 1].parent_i = i;
                cellDetails[i - 1][j + 1].parent_j = j;
                printf("The destination cell is found\n");
                tracePath(cellDetails, dest_row, dest_col, out_grid);
                foundDest = true;
                return;
            }
            else if (closedList[i - 1][j + 1] == false && isUnBlocked(grid, i - 1, j + 1) == true)
            {
                gNew = cellDetails[i][j].g + 1.414;
                hNew = calculateHValue(i - 1, j + 1, dest_row, dest_col);
                fNew = gNew + hNew;
                if (cellDetails[i - 1][j + 1].f == FLT_MAX || cellDetails[i - 1][j + 1].f > fNew)
                {
                    openListSize += 1;
                    openList[openListSize - 1] = newOpenCellData(i - 1, j + 1, fNew);
                    cellDetails[i - 1][j + 1] = newCell(i, j, fNew, gNew, hNew);
                }
            }
        }
 
        //----------- 6th Successor (North-West) ------------
        if (isValid(i - 1, j - 1) == true)
        {
            if (isDestination(i - 1, j - 1, dest_row, dest_col) == true)
            {
                cellDetails[i - 1][j - 1].parent_i = i;
                cellDetails[i - 1][j - 1].parent_j = j;
                printf("The destination cell is found\n");
                tracePath(cellDetails, dest_row, dest_col, out_grid);
                foundDest = true;
                return;
            }
            else if (closedList[i - 1][j - 1] == false && isUnBlocked(grid, i - 1, j - 1) == true)
            {
                gNew = cellDetails[i][j].g + 1.414;
                hNew = calculateHValue(i - 1, j - 1, dest_row, dest_col);
                fNew = gNew + hNew;
                if (cellDetails[i - 1][j - 1].f == FLT_MAX || cellDetails[i - 1][j - 1].f > fNew)
                {
                    openListSize += 1;
                    openList[openListSize - 1] = newOpenCellData(i - 1, j - 1, fNew);
                    cellDetails[i - 1][j - 1] = newCell(i, j, fNew, gNew, hNew);
                }
            }
        }
 
        //----------- 7th Successor (South-East) ------------
        if (isValid(i + 1, j + 1) == true)
        {
            if (isDestination(i + 1, j + 1, dest_row, dest_col) == true)
            {
                cellDetails[i + 1][j + 1].parent_i = i;
                cellDetails[i + 1][j + 1].parent_j = j;
                printf("The destination cell is found\n");
                tracePath(cellDetails, dest_row, dest_col, out_grid);
                foundDest = true;
                return;
            }
            else if (closedList[i + 1][j + 1] == false && isUnBlocked(grid, i + 1, j + 1) == true) 
            {
                gNew = cellDetails[i][j].g + 1.414;
                hNew = calculateHValue(i + 1, j + 1, dest_row, dest_col);
                fNew = gNew + hNew;
                if (cellDetails[i + 1][j + 1].f == FLT_MAX || cellDetails[i + 1][j + 1].f > fNew)
                {
                    openListSize += 1;
                    openList[openListSize - 1] = newOpenCellData(i + 1, j + 1, fNew);
                    cellDetails[i + 1][j + 1] = newCell(i, j, fNew, gNew, hNew);
                }
            }
        }
 
        //----------- 8th Successor (South-West) ------------
        if (isValid(i + 1, j - 1) == true)
        {
            if (isDestination(i + 1, j - 1, dest_row, dest_col) == true)
            {
                cellDetails[i + 1][j - 1].parent_i = i;
                cellDetails[i + 1][j - 1].parent_j = j;
                printf("The destination cell is found\n");
                tracePath(cellDetails, dest_row, dest_col, out_grid);
                foundDest = true;
                return;
            }
            else if (closedList[i + 1][j - 1] == false && isUnBlocked(grid, i + 1, j - 1) == true)
            {
                gNew = cellDetails[i][j].g + 1.414;
                hNew = calculateHValue(i + 1, j - 1, dest_row, dest_col);
                fNew = gNew + hNew;
                if (cellDetails[i + 1][j - 1].f == FLT_MAX || cellDetails[i + 1][j - 1].f > fNew)
                {
                    openListSize += 1;
                    openList[openListSize - 1] = newOpenCellData(i + 1, j - 1, fNew);
                    cellDetails[i + 1][j - 1] = newCell(i, j, fNew, gNew, hNew);
                }
            }
        }
    }
    if (foundDest == false)
        printf("Failed to find the Destination Cell\n");
 
    return;
}
 
int main(int argc, char** argv) {
     
    int rank, size;
    double tstart, tend;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int row = ROW;
    int col = COL;
    int chunk = ceil((double)NUM_SEARCH / size);
    int start = rank * chunk;
    int end = (rank+1) * chunk - 1;
    if (end > NUM_SEARCH)
        end = NUM_SEARCH - 1;

    int** grid = allocate_matrix(row, col);
    int** out_grid = allocate_matrix(row, col);
    if (rank == 0)
    {
        printf("Reading file\n");
        readmatrix(row, col, grid, "in_path.txt");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < ROW; i++)
        for (int j = 0; j < COL; j++)
            MPI_Bcast(&grid[i][j], 1, MPI_INT, 0, MPI_COMM_WORLD);

    tstart = MPI_Wtime();
    

    if(start <= NUM_SEARCH)
    {
        for (int i = start; i <= end; i++)
        {
            for (int j = 0; j < ROW; j++)
                for (int k = 0; k < COL; k++)
                    out_grid[j][k] = grid[j][k];

            printf("Start %d a* searching in %d thread\n", i, rank);
            aStar(grid, out_grid, 0, 0, ROW-1, COL-1);

            char filename[100];
            sprintf(filename, "out_path%d.txt", i);
            writematrix(ROW, COL, out_grid, filename);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    tend = MPI_Wtime();
    if (rank == 0)
    {
        double time_taken_parallel = tend - tstart;
        printf("Time spend = %f\n", time_taken_parallel);
    }
    MPI_Finalize();
    return 0;
}