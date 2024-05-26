// #include <mpi.h>
#include <stdio.h>
#include </home/daniil/Parallel/vec.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#define ROW 5
#define COL 5

typedef struct {
    int first, second;
}Pair;

Pair newPair(int first, int second) {
    Pair p;
    p.first = first;
    p.second = second;
    return p;
}

void push(Pair* array, int* size, Pair newData) 
{
    *size += 1;
    Pair* tmp = (Pair*)realloc(array, *size * sizeof(Pair));
    array = (Pair*)realloc(array, (*size) * sizeof(Pair));

    for (int i = 1; i < *size; i++)
    {
        array[i] = array[i + 1];
    }
    array[0] = newData;
    return;
}

// void push(Pair* array, int* size, Pair newData) 
// {
//     // Увеличиваем размер массива на 1
//     Pair* tmp = (Pair*)malloc(*size * sizeof(Pair));
//     for(int i = 0; i < *size; i++)
//     {        
//         tmp[i] = array[i];
//     }

//     *size += 1;
//     array = (Pair*)realloc(array, (*size) * sizeof(Pair));

//     for(int i = 1; i < *size; i++)
//     {        
//         array[i] = tmp[i - 1];
//     }
//     free(tmp);
//     // printf("push size %d\n", *size);
//     array[0] = newData;
//     for(int i = 0; i < *size; i++)
//     {        
//         printf("%d %d\n", array[i].first, array[i].second);
//     }
//     printf("was in push\n");
//     return;
// }

Pair pop(Pair* array, int* size) 
{
    Pair popValue = array[0];
    *size -= 1; 
    for(int i = 0; i < *size - 1; i++)
            array[i] = array[i + 1];

    if (*size >= 0)
    {
        Pair* tmp = (Pair*)realloc(array, *size * sizeof(Pair));
        array = tmp;
        tmp = NULL;
    }

    return popValue;
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

void append(OpenCellData* array, int* size, OpenCellData newData) {
    // Увеличиваем размер массива на 1
    *size += 1;
    array = (OpenCellData*)realloc(array, (*size) * sizeof(OpenCellData));
    array[*size - 1] = newData;
    return;
}

OpenCellData* remove_element(OpenCellData* array, int index, int* size){
    int i;
    for(i = index; i < *size - 1; i++)
        array[i] = array[i + 1];
    array = realloc(array, (*size - 1) * sizeof(OpenCellData));
    *size -= 1;
    return array;
}

typedef struct{
    // Row and Column index of its parent
    int parent_i, parent_j;
    // f = g + h
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

void tracePath(cell** cellDetails, int dest_row, int dest_col)
{
    printf("\nThe Path is \n");
    int row = dest_row;
    int col = dest_col;
    Pair* Path = (Pair *)malloc(0);
    int size = 0;
    Pair tmpP;
    while (!(cellDetails[row][col].parent_i == row && cellDetails[row][col].parent_j == col)) 
    {
        tmpP = newPair(row, col);
        printf("%d %d\n", tmpP.first, tmpP.second);
        size += 1;
        Path = (Pair*)realloc(Path, size * sizeof(Pair));

        for (int i = size - 1; i > 0; i--)
        {
            Path[i] = Path[i - 1];
        }
        Path[0] = tmpP;
        
        int temp_row = cellDetails[row][col].parent_i;
        int temp_col = cellDetails[row][col].parent_j;
        row = temp_row;
        col = temp_col;
        if (cellDetails[row][col].parent_i == row && cellDetails[row][col].parent_j == col)
        {
            printf("\n");
            for (int i = 0; i < size; i++)
            {
                printf("%d %d\n", Path[i].first, Path[i].second);
            }
            tmpP = newPair(row, col);
            size += 1;
            Path = (Pair*)realloc(Path, size * sizeof(Pair));

            for (int i = size - 1; i > 0; i--)
            {
                Path[i] = Path[i - 1];
            }
            Path[0] = tmpP;

            printf("\n");
            for (int i = 0; i < size; i++)
            {
                printf("%d %d\n", Path[i].first, Path[i].second);
            }
        }
    }

    // tmpP = newPair(row, col);
    // printf("%d %d\n", tmpP.first, tmpP.second);
    // size += 1;
    // Path = (Pair*)realloc(Path, size * sizeof(Pair));
    // printf("\n");
    // printf("\n");
    // for (int i = size - 1; i > 0; i--)
    // {
    //     Path[i] = Path[i - 1];
    //     printf("%d %d\n", Path[i].first, Path[i].second);
    // }
    // Path[0] = tmpP;
    printf("\n");
    printf("\n");
    for (int i = 0; i < size; i++)
    {
        printf("%d %d\n", Path[i].first, Path[i].second);
    }
    printf("\n");
    while (size > 0) {
        Pair p = Path[0];
        size -= 1; 
        for(int i = 0; i < size; i++)
            Path[i] = Path[i + 1];
        Path = (Pair*)realloc(Path, size * sizeof(Pair));
        printf("-> (%d,%d) \n", p.first, p.second);
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

void aTest(int** grid, int src_row, int src_col, int dest_row, int dest_col)
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
    OpenCellData* openList = (OpenCellData *)malloc(sizeof(OpenCellData));
    openList[0] = newOpenCellData(i, j, 0.0);
    int openListSize = 1;
    // We set this boolean value as false as initially
    // the destination is not reached.
    bool foundDest = false;

    OpenCellData tmpOpenCell;

    while (openListSize > 0) {
        OpenCellData p = openList[0];
 
        // Remove this vertex from the open list
        remove_element(openList, 0, &openListSize);
 
        // Add this vertex to the closed list
        i = p.row;
        j = p.col;
        closedList[i][j] = true;
 
        /*
         Generating all the 8 successor of this cell
 
             N.W   N   N.E
               \   |   /
                \  |  /
             W----Cell----E
                  / | \
                /   |  \
             S.W    S   S.E
 
         Cell-->Popped Cell (i, j)
         N -->  North       (i-1, j)
         S -->  South       (i+1, j)
         E -->  East        (i, j+1)
         W -->  West           (i, j-1)
         N.E--> North-East  (i-1, j+1)
         N.W--> North-West  (i-1, j-1)
         S.E--> South-East  (i+1, j+1)
         S.W--> South-West  (i+1, j-1)*/
 
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
                tracePath(cellDetails, dest_row, dest_col);
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
                    tmpOpenCell = newOpenCellData(i - 1, j, fNew);
                    append(openList, &openListSize, tmpOpenCell);
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
                tracePath(cellDetails, dest_row, dest_col);
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
                    tmpOpenCell = newOpenCellData(i + 1, j, fNew);
                    append(openList, &openListSize, tmpOpenCell);
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
                tracePath(cellDetails, dest_row, dest_col);
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
                    tmpOpenCell = newOpenCellData(i, j + 1, fNew);
                    append(openList, &openListSize, tmpOpenCell);
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
                tracePath(cellDetails, dest_row, dest_col);
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
                    tmpOpenCell = newOpenCellData(i, j - 1, fNew);
                    append(openList, &openListSize, tmpOpenCell);
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
                tracePath(cellDetails, dest_row, dest_col);
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
                    tmpOpenCell = newOpenCellData(i - 1, j + 1, fNew);
                    append(openList, &openListSize, tmpOpenCell);
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
                tracePath(cellDetails, dest_row, dest_col);
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
                    tmpOpenCell = newOpenCellData(i - 1, j - 1, fNew);
                    append(openList, &openListSize, tmpOpenCell);
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
                tracePath(cellDetails, dest_row, dest_col);
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

                    tmpOpenCell = newOpenCellData(i + 1, j + 1, fNew);
                    append(openList, &openListSize, tmpOpenCell);
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
                tracePath(cellDetails, dest_row, dest_col);
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
                    tmpOpenCell = newOpenCellData(i + 1, j - 1, fNew);
                    append(openList, &openListSize, tmpOpenCell);
                    cellDetails[i + 1][j - 1] = newCell(i, j, fNew, gNew, hNew);
                }
            }
        }
    }

    //     // Only process this cell if this is a valid one
    //     if (isValid(i + 1, j - 1) == true) {
    //         // If the destination cell is the same as the
    //         // current successor
    //         if (isDestination(i + 1, j - 1, dest) == true) {
    //             // Set the Parent of the destination cell
    //             cellDetails[i + 1][j - 1].parent_i = i;
    //             cellDetails[i + 1][j - 1].parent_j = j;
    //             printf("The destination cell is found\n");
    //             tracePath(cellDetails, dest);
    //             foundDest = true;
    //             return;
    //         }
 
    //         // If the successor is already on the closed
    //         // list or if it is blocked, then ignore it.
    //         // Else do the following
    //         else if (closedList[i + 1][j - 1] == false
    //                  && isUnBlocked(grid, i + 1, j - 1)
    //                         == true) {
    //             gNew = cellDetails[i][j].g + 1.414;
    //             hNew = calculateHValue(i + 1, j - 1, dest);
    //             fNew = gNew + hNew;
 
    //             // If it isn’t on the open list, add it to
    //             // the open list. Make the current square
    //             // the parent of this square. Record the
    //             // f, g, and h costs of the square cell
    //             //                OR
    //             // If it is on the open list already, check
    //             // to see if this path to that square is
    //             // better, using 'f' cost as the measure.
    //             if (cellDetails[i + 1][j - 1].f == FLT_MAX
    //                 || cellDetails[i + 1][j - 1].f > fNew) {
    //                 openList.insert(make_pair(
    //                     fNew, make_pair(i + 1, j - 1)));
 
    //                 // Update the details of this cell
    //                 cellDetails[i + 1][j - 1].f = fNew;
    //                 cellDetails[i + 1][j - 1].g = gNew;
    //                 cellDetails[i + 1][j - 1].h = hNew;
    //                 cellDetails[i + 1][j - 1].parent_i = i;
    //                 cellDetails[i + 1][j - 1].parent_j = j;
    //             }
    //         }
    //     }
    // }
 
    if (foundDest == false)
        printf("Failed to find the Destination Cell\n");
 
    return;
}



















// void aStarSearch(int** grid, int src_row, int src_col, int dest_row, int dest_col)
// {
//     // If the source is out of range
//     if (isValid(src_row, src_col) == false) {
//         printf("Source is invalid\n");
//         return;
//     }
 
//     // If the destination is out of range
//     if (isValid(dest_row, dest_col) == false) {
//         printf("Destination is invalid\n");
//         return;
//     }
 
//     // Either the source or the destination is blocked
//     if (isUnBlocked(grid, src_row, src_col) == false
//         || isUnBlocked(grid, dest_row, dest_col)
//                == false) {
//         printf("Source or the destination is blocked\n");
//         return;
//     }
 
//     // If the destination cell is the same as source cell
//     if (isDestination(src_row, src_col, dest)
//         == true) {
//         printf("We are already at the destination\n");
//         return;
//     }
 
//     // Create a closed list and initialise it to false which
//     // means that no cell has been included yet This closed
//     // list is implemented as a boolean 2D array

//     bool closedList[ROW][COL];
//     memset(closedList, false, sizeof(closedList));
 
//     // Declare a 2D array of structure to hold the details
//     // of that cell
//     struct cell** cellDetails[ROW][COL];
 
//     int i, j;
 
//     for (i = 0; i < ROW; i++) {
//         for (j = 0; j < COL; j++) {
//             cellDetails[i][j].f = FLT_MAX;
//             cellDetails[i][j].g = FLT_MAX;
//             cellDetails[i][j].h = FLT_MAX;
//             cellDetails[i][j].parent_i = -1;
//             cellDetails[i][j].parent_j = -1;
//         }
//     }
 
//     // Initialising the parameters of the starting node
//     i = src_row, j = src_col;
//     cellDetails[i][j].f = 0.0;
//     cellDetails[i][j].g = 0.0;
//     cellDetails[i][j].h = 0.0;
//     cellDetails[i][j].parent_i = i;
//     cellDetails[i][j].parent_j = j;
 
//     /*
//      Create an open list having information as-
//      <f, <i, j>>
//      where f = g + h,
//      and i, j are the row and column index of that cell
//      Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
//      This open list is implemented as a set of pair of
//      pair.*/
//     struct openCellData* openList[];
 
//     // Put the starting cell on the open list and set its
//     // 'f' as 0+
//     openList.insert(make_pair(0.0, make_pair(i, j)));
 
//     // We set this boolean value as false as initially
//     // the destination is not reached.
//     bool foundDest = false;
 
//     while (!openList.empty()) {
//         openCellData p = *openList.begin();
 
//         // Remove this vertex from the open list
//         openList.erase(openList.begin());
 
//         // Add this vertex to the closed list
//         i = p.second.first;
//         j = p.second.second;
//         closedList[i][j] = true;
 
//         /*
//          Generating all the 8 successor of this cell
 
//              N.W   N   N.E
//                \   |   /
//                 \  |  /
//              W----Cell----E
//                   / | \
//                 /   |  \
//              S.W    S   S.E
 
//          Cell-->Popped Cell (i, j)
//          N -->  North       (i-1, j)
//          S -->  South       (i+1, j)
//          E -->  East        (i, j+1)
//          W -->  West           (i, j-1)
//          N.E--> North-East  (i-1, j+1)
//          N.W--> North-West  (i-1, j-1)
//          S.E--> South-East  (i+1, j+1)
//          S.W--> South-West  (i+1, j-1)*/
 
//         // To store the 'g', 'h' and 'f' of the 8 successors
//         double gNew, hNew, fNew;
 
//         //----------- 1st Successor (North) ------------
 
//         // Only process this cell if this is a valid one
//         if (isValid(i - 1, j) == true) {
//             // If the destination cell is the same as the
//             // current successor
//             if (isDestination(i - 1, j, dest) == true) {
//                 // Set the Parent of the destination cell
//                 cellDetails[i - 1][j].parent_i = i;
//                 cellDetails[i - 1][j].parent_j = j;
//                 printf("The destination cell is found\n");
//                 tracePath(cellDetails, dest);
//                 foundDest = true;
//                 return;
//             }
//             // If the successor is already on the closed
//             // list or if it is blocked, then ignore it.
//             // Else do the following
//             else if (closedList[i - 1][j] == false
//                      && isUnBlocked(grid, i - 1, j)
//                             == true) {
//                 gNew = cellDetails[i][j].g + 1.0;
//                 hNew = calculateHValue(i - 1, j, dest);
//                 fNew = gNew + hNew;
 
//                 // If it isn’t on the open list, add it to
//                 // the open list. Make the current square
//                 // the parent of this square. Record the
//                 // f, g, and h costs of the square cell
//                 //                OR
//                 // If it is on the open list already, check
//                 // to see if this path to that square is
//                 // better, using 'f' cost as the measure.
//                 if (cellDetails[i - 1][j].f == FLT_MAX
//                     || cellDetails[i - 1][j].f > fNew) {
//                     openList.insert(make_pair(
//                         fNew, make_pair(i - 1, j)));
 
//                     // Update the details of this cell
//                     cellDetails[i - 1][j].f = fNew;
//                     cellDetails[i - 1][j].g = gNew;
//                     cellDetails[i - 1][j].h = hNew;
//                     cellDetails[i - 1][j].parent_i = i;
//                     cellDetails[i - 1][j].parent_j = j;
//                 }
//             }
//         }
 
//         //----------- 2nd Successor (South) ------------
 
//         // Only process this cell if this is a valid one
//         if (isValid(i + 1, j) == true) {
//             // If the destination cell is the same as the
//             // current successor
//             if (isDestination(i + 1, j, dest) == true) {
//                 // Set the Parent of the destination cell
//                 cellDetails[i + 1][j].parent_i = i;
//                 cellDetails[i + 1][j].parent_j = j;
//                 printf("The destination cell is found\n");
//                 tracePath(cellDetails, dest);
//                 foundDest = true;
//                 return;
//             }
//             // If the successor is already on the closed
//             // list or if it is blocked, then ignore it.
//             // Else do the following
//             else if (closedList[i + 1][j] == false
//                      && isUnBlocked(grid, i + 1, j)
//                             == true) {
//                 gNew = cellDetails[i][j].g + 1.0;
//                 hNew = calculateHValue(i + 1, j, dest);
//                 fNew = gNew + hNew;
 
//                 // If it isn’t on the open list, add it to
//                 // the open list. Make the current square
//                 // the parent of this square. Record the
//                 // f, g, and h costs of the square cell
//                 //                OR
//                 // If it is on the open list already, check
//                 // to see if this path to that square is
//                 // better, using 'f' cost as the measure.
//                 if (cellDetails[i + 1][j].f == FLT_MAX
//                     || cellDetails[i + 1][j].f > fNew) {
//                     openList.insert(make_pair(
//                         fNew, make_pair(i + 1, j)));
//                     // Update the details of this cell
//                     cellDetails[i + 1][j].f = fNew;
//                     cellDetails[i + 1][j].g = gNew;
//                     cellDetails[i + 1][j].h = hNew;
//                     cellDetails[i + 1][j].parent_i = i;
//                     cellDetails[i + 1][j].parent_j = j;
//                 }
//             }
//         }
 
//         //----------- 3rd Successor (East) ------------
 
//         // Only process this cell if this is a valid one
//         if (isValid(i, j + 1) == true) {
//             // If the destination cell is the same as the
//             // current successor
//             if (isDestination(i, j + 1, dest) == true) {
//                 // Set the Parent of the destination cell
//                 cellDetails[i][j + 1].parent_i = i;
//                 cellDetails[i][j + 1].parent_j = j;
//                 printf("The destination cell is found\n");
//                 tracePath(cellDetails, dest);
//                 foundDest = true;
//                 return;
//             }
 
//             // If the successor is already on the closed
//             // list or if it is blocked, then ignore it.
//             // Else do the following
//             else if (closedList[i][j + 1] == false
//                      && isUnBlocked(grid, i, j + 1)
//                             == true) {
//                 gNew = cellDetails[i][j].g + 1.0;
//                 hNew = calculateHValue(i, j + 1, dest);
//                 fNew = gNew + hNew;
 
//                 // If it isn’t on the open list, add it to
//                 // the open list. Make the current square
//                 // the parent of this square. Record the
//                 // f, g, and h costs of the square cell
//                 //                OR
//                 // If it is on the open list already, check
//                 // to see if this path to that square is
//                 // better, using 'f' cost as the measure.
//                 if (cellDetails[i][j + 1].f == FLT_MAX
//                     || cellDetails[i][j + 1].f > fNew) {
//                     openList.insert(make_pair(
//                         fNew, make_pair(i, j + 1)));
 
//                     // Update the details of this cell
//                     cellDetails[i][j + 1].f = fNew;
//                     cellDetails[i][j + 1].g = gNew;
//                     cellDetails[i][j + 1].h = hNew;
//                     cellDetails[i][j + 1].parent_i = i;
//                     cellDetails[i][j + 1].parent_j = j;
//                 }
//             }
//         }
 
//         //----------- 4th Successor (West) ------------
 
//         // Only process this cell if this is a valid one
//         if (isValid(i, j - 1) == true) {
//             // If the destination cell is the same as the
//             // current successor
//             if (isDestination(i, j - 1, dest) == true) {
//                 // Set the Parent of the destination cell
//                 cellDetails[i][j - 1].parent_i = i;
//                 cellDetails[i][j - 1].parent_j = j;
//                 printf("The destination cell is found\n");
//                 tracePath(cellDetails, dest);
//                 foundDest = true;
//                 return;
//             }
 
//             // If the successor is already on the closed
//             // list or if it is blocked, then ignore it.
//             // Else do the following
//             else if (closedList[i][j - 1] == false
//                      && isUnBlocked(grid, i, j - 1)
//                             == true) {
//                 gNew = cellDetails[i][j].g + 1.0;
//                 hNew = calculateHValue(i, j - 1, dest);
//                 fNew = gNew + hNew;
 
//                 // If it isn’t on the open list, add it to
//                 // the open list. Make the current square
//                 // the parent of this square. Record the
//                 // f, g, and h costs of the square cell
//                 //                OR
//                 // If it is on the open list already, check
//                 // to see if this path to that square is
//                 // better, using 'f' cost as the measure.
//                 if (cellDetails[i][j - 1].f == FLT_MAX
//                     || cellDetails[i][j - 1].f > fNew) {
//                     openList.insert(make_pair(
//                         fNew, make_pair(i, j - 1)));
 
//                     // Update the details of this cell
//                     cellDetails[i][j - 1].f = fNew;
//                     cellDetails[i][j - 1].g = gNew;
//                     cellDetails[i][j - 1].h = hNew;
//                     cellDetails[i][j - 1].parent_i = i;
//                     cellDetails[i][j - 1].parent_j = j;
//                 }
//             }
//         }
 
//         //----------- 5th Successor (North-East)
//         //------------
 
//         // Only process this cell if this is a valid one
//         if (isValid(i - 1, j + 1) == true) {
//             // If the destination cell is the same as the
//             // current successor
//             if (isDestination(i - 1, j + 1, dest) == true) {
//                 // Set the Parent of the destination cell
//                 cellDetails[i - 1][j + 1].parent_i = i;
//                 cellDetails[i - 1][j + 1].parent_j = j;
//                 printf("The destination cell is found\n");
//                 tracePath(cellDetails, dest);
//                 foundDest = true;
//                 return;
//             }
 
//             // If the successor is already on the closed
//             // list or if it is blocked, then ignore it.
//             // Else do the following
//             else if (closedList[i - 1][j + 1] == false
//                      && isUnBlocked(grid, i - 1, j + 1)
//                             == true) {
//                 gNew = cellDetails[i][j].g + 1.414;
//                 hNew = calculateHValue(i - 1, j + 1, dest);
//                 fNew = gNew + hNew;
 
//                 // If it isn’t on the open list, add it to
//                 // the open list. Make the current square
//                 // the parent of this square. Record the
//                 // f, g, and h costs of the square cell
//                 //                OR
//                 // If it is on the open list already, check
//                 // to see if this path to that square is
//                 // better, using 'f' cost as the measure.
//                 if (cellDetails[i - 1][j + 1].f == FLT_MAX
//                     || cellDetails[i - 1][j + 1].f > fNew) {
//                     openList.insert(make_pair(
//                         fNew, make_pair(i - 1, j + 1)));
 
//                     // Update the details of this cell
//                     cellDetails[i - 1][j + 1].f = fNew;
//                     cellDetails[i - 1][j + 1].g = gNew;
//                     cellDetails[i - 1][j + 1].h = hNew;
//                     cellDetails[i - 1][j + 1].parent_i = i;
//                     cellDetails[i - 1][j + 1].parent_j = j;
//                 }
//             }
//         }
 
//         //----------- 6th Successor (North-West)
//         //------------
 
//         // Only process this cell if this is a valid one
//         if (isValid(i - 1, j - 1) == true) {
//             // If the destination cell is the same as the
//             // current successor
//             if (isDestination(i - 1, j - 1, dest) == true) {
//                 // Set the Parent of the destination cell
//                 cellDetails[i - 1][j - 1].parent_i = i;
//                 cellDetails[i - 1][j - 1].parent_j = j;
//                 printf("The destination cell is found\n");
//                 tracePath(cellDetails, dest);
//                 foundDest = true;
//                 return;
//             }
 
//             // If the successor is already on the closed
//             // list or if it is blocked, then ignore it.
//             // Else do the following
//             else if (closedList[i - 1][j - 1] == false
//                      && isUnBlocked(grid, i - 1, j - 1)
//                             == true) {
//                 gNew = cellDetails[i][j].g + 1.414;
//                 hNew = calculateHValue(i - 1, j - 1, dest);
//                 fNew = gNew + hNew;
 
//                 // If it isn’t on the open list, add it to
//                 // the open list. Make the current square
//                 // the parent of this square. Record the
//                 // f, g, and h costs of the square cell
//                 //                OR
//                 // If it is on the open list already, check
//                 // to see if this path to that square is
//                 // better, using 'f' cost as the measure.
//                 if (cellDetails[i - 1][j - 1].f == FLT_MAX
//                     || cellDetails[i - 1][j - 1].f > fNew) {
//                     openList.insert(make_pair(
//                         fNew, make_pair(i - 1, j - 1)));
//                     // Update the details of this cell
//                     cellDetails[i - 1][j - 1].f = fNew;
//                     cellDetails[i - 1][j - 1].g = gNew;
//                     cellDetails[i - 1][j - 1].h = hNew;
//                     cellDetails[i - 1][j - 1].parent_i = i;
//                     cellDetails[i - 1][j - 1].parent_j = j;
//                 }
//             }
//         }
 
//         //----------- 7th Successor (South-East)
//         //------------
 
//         // Only process this cell if this is a valid one
//         if (isValid(i + 1, j + 1) == true) {
//             // If the destination cell is the same as the
//             // current successor
//             if (isDestination(i + 1, j + 1, dest) == true) {
//                 // Set the Parent of the destination cell
//                 cellDetails[i + 1][j + 1].parent_i = i;
//                 cellDetails[i + 1][j + 1].parent_j = j;
//                 printf("The destination cell is found\n");
//                 tracePath(cellDetails, dest);
//                 foundDest = true;
//                 return;
//             }
 
//             // If the successor is already on the closed
//             // list or if it is blocked, then ignore it.
//             // Else do the following
//             else if (closedList[i + 1][j + 1] == false
//                      && isUnBlocked(grid, i + 1, j + 1)
//                             == true) {
//                 gNew = cellDetails[i][j].g + 1.414;
//                 hNew = calculateHValue(i + 1, j + 1, dest);
//                 fNew = gNew + hNew;
 
//                 // If it isn’t on the open list, add it to
//                 // the open list. Make the current square
//                 // the parent of this square. Record the
//                 // f, g, and h costs of the square cell
//                 //                OR
//                 // If it is on the open list already, check
//                 // to see if this path to that square is
//                 // better, using 'f' cost as the measure.
//                 if (cellDetails[i + 1][j + 1].f == FLT_MAX
//                     || cellDetails[i + 1][j + 1].f > fNew) {
//                     openList.insert(make_pair(
//                         fNew, make_pair(i + 1, j + 1)));
 
//                     // Update the details of this cell
//                     cellDetails[i + 1][j + 1].f = fNew;
//                     cellDetails[i + 1][j + 1].g = gNew;
//                     cellDetails[i + 1][j + 1].h = hNew;
//                     cellDetails[i + 1][j + 1].parent_i = i;
//                     cellDetails[i + 1][j + 1].parent_j = j;
//                 }
//             }
//         }
 
//         //----------- 8th Successor (South-West)
//         //------------
 
//         // Only process this cell if this is a valid one
//         if (isValid(i + 1, j - 1) == true) {
//             // If the destination cell is the same as the
//             // current successor
//             if (isDestination(i + 1, j - 1, dest) == true) {
//                 // Set the Parent of the destination cell
//                 cellDetails[i + 1][j - 1].parent_i = i;
//                 cellDetails[i + 1][j - 1].parent_j = j;
//                 printf("The destination cell is found\n");
//                 tracePath(cellDetails, dest);
//                 foundDest = true;
//                 return;
//             }
 
//             // If the successor is already on the closed
//             // list or if it is blocked, then ignore it.
//             // Else do the following
//             else if (closedList[i + 1][j - 1] == false
//                      && isUnBlocked(grid, i + 1, j - 1)
//                             == true) {
//                 gNew = cellDetails[i][j].g + 1.414;
//                 hNew = calculateHValue(i + 1, j - 1, dest);
//                 fNew = gNew + hNew;
 
//                 // If it isn’t on the open list, add it to
//                 // the open list. Make the current square
//                 // the parent of this square. Record the
//                 // f, g, and h costs of the square cell
//                 //                OR
//                 // If it is on the open list already, check
//                 // to see if this path to that square is
//                 // better, using 'f' cost as the measure.
//                 if (cellDetails[i + 1][j - 1].f == FLT_MAX
//                     || cellDetails[i + 1][j - 1].f > fNew) {
//                     openList.insert(make_pair(
//                         fNew, make_pair(i + 1, j - 1)));
 
//                     // Update the details of this cell
//                     cellDetails[i + 1][j - 1].f = fNew;
//                     cellDetails[i + 1][j - 1].g = gNew;
//                     cellDetails[i + 1][j - 1].h = hNew;
//                     cellDetails[i + 1][j - 1].parent_i = i;
//                     cellDetails[i + 1][j - 1].parent_j = j;
//                 }
//             }
//         }
//     }
 
//     // When the destination cell is not found and the open
//     // list is empty, then we conclude that we failed to
//     // reach the destination cell. This may happen when the
//     // there is no way to destination cell (due to
//     // blockages)
//     if (foundDest == false)
//         printf("Failed to find the Destination Cell\n");
 
//     return;
// }
 
int main(int argc, char** argv) {
    int rank, size;
    // MPI_Init(&argc, &argv);
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &size);
    int** grid = allocate_matrix(5, 5);

    readmatrix(ROW, COL, grid, "in_path.txt");

    // int* ar = (int*)malloc(5 * sizeof(int));
    // for (int i = 0; i < 5; i++)
    //     ar[i] = i;

    // for (int i = 0; i < 5; i++)
    //     printf("%d ", ar[i]);
    // printf("\n");

    // ar = (int*)realloc(ar, 10 * sizeof(int));
    // for (int i = 0; i < 10; i++)
    //     printf("%d ", ar[i]);
    // printf("\n");
    // for (int i = 0; i < 5; i++)
    // {
    //     for (int j = 0; j < 5; j++)
    //     {
    //         printf("%d", grid[i][j]);
    //     }
    //     printf("\n");
    // }

    aTest(grid, 0, 0, 4, 4);
    // MPI_Finalize();
    return 0;
}