#include "combinatorics.h"
#include "FpVector.h"

typedef struct {
    uint p;
    uint rows;
    uint columns;
    Vector **vectors;
} Matrix;

uint Matrix_getSize(uint p, uint rows, uint cols);
Matrix *Matrix_initialize(char *memory, uint p, uint rows, uint cols);
Matrix *Matrix_construct(uint p, uint rows, uint cols);
void Matrix_free(Matrix *M);

uint Matrix_toString(char *buffer, Matrix *M);
void Matrix_print(Matrix *matrix);
uint Matrix_getSliceSize(uint p, uint rows);
Matrix *Matrix_slice(Matrix *M, char *memory, uint row_min, uint row_max, uint column_min, uint column_max);
void Matrix_printSlice(Matrix *M, uint col_end, uint col_start);

void Matrix_serialize(char **buffer, Matrix *v);
Matrix *Matrix_deserialize(char **buffer);

void Matrix_getRowPermutation(Matrix *M, uint *result);
void Matrix_applyRowPermutation(Matrix *M, uint *permutation, uint rows);

// Row reduce M. For each column i, column_to_pivot_row[i] is equal to: the row
// with a pivot in column i or -1 if no such row exists.
// When you delete the -1's from column_to_pivot_row, it looks like [0,1,2,...,rank(M)].
void rowReduce(Matrix *M, int *column_to_pivot_row, uint, uint);

typedef struct {
    Matrix *matrix;
    int *column_to_pivot_row;
} Subspace;

size_t Subspace_getSize(uint p, uint rows, uint columns);

Subspace *Subspace_construct(uint p, uint rows, uint columns);
void Subspace_free(Subspace *k);

void Subspace_serialize(char **buffer, Subspace *subspace);
Subspace *Subspace_deserialize(char **buffer);

Subspace *Matrix_computeKernel(Matrix *M, int *column_to_pivot_row, uint first_source_column);
uint Matrix_extendImage(Matrix *matrix, uint first_source_column, uint source_dimension, uint first_empty_row, int *current_pivots, Subspace *desired_image);

void Matrix_quasiInverse_apply(Vector *target, Subspace *image, Matrix *quasi_inverse, Vector *input);