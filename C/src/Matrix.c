#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "Matrix.h"


uint Matrix_getSize(uint p, uint rows, uint cols){
    assert(cols < MAX_DIMENSION);
    return sizeof(Matrix) 
      + rows * (sizeof(Vector*) + Vector_getSize(p, cols, 0));
}

Matrix *Matrix_initialize(char *memory, uint p, uint rows, uint columns)  {
    Matrix *matrix = (Matrix*)memory;
    Vector **vector_ptr = (Vector**)(matrix+1);
    char *container_ptr = (char*)(vector_ptr + rows); //
    matrix->p = p;
    matrix->rows = rows;
    matrix->columns = columns;
    matrix->vectors = vector_ptr;
    for(uint row = 0; row < rows; row++){
        *vector_ptr = Vector_initialize(p, &container_ptr, columns, 0);
        vector_ptr ++; 
    }
    return matrix;
}

void Matrix_serialize(char **buffer, Matrix *M){
    size_t size = Matrix_getSize(M->p, M->rows, M->columns);
    memcpy(*buffer, M, sizeof(Matrix));
    *buffer += sizeof(Matrix);
    *buffer += M->rows * sizeof(Vector*);
    for(uint row = 0; row < M->rows; row++){
        Vector_serialize(buffer, M->vectors[row]);
    }
}

Matrix *Matrix_deserialize(char **buffer){
    Matrix *M = (Matrix*)*buffer;
    char *start_ptr = *buffer;
    *buffer += sizeof(Matrix);
    Vector **vector_ptr = (Vector**)*buffer;
    M->vectors = vector_ptr;
    uint rows = M->rows;    
    *buffer += rows * sizeof(Vector*);
    for(uint row = 0; row < rows; row++){
        *vector_ptr = Vector_deserialize(M->p, buffer);
        vector_ptr ++;
    }
    assert(start_ptr + Matrix_getSize(M->p, M->rows, M->columns) == *buffer);
    return M;
}

Matrix *Matrix_construct(uint p,  uint rows, uint columns)  {
    char *M = malloc(Matrix_getSize(p, rows, columns));
    // printf("columns: %d, rows: %d\n", columns, rows);
    return Matrix_initialize(M, p, rows, columns);
}

void Matrix_free(Matrix *M){
    free(M);
}

Vector *Matrix_getRow(Matrix *M, uint row){
    assert(row < M->rows);
    return M->vectors[row];
}


uint Matrix_getSliceSize(uint p __attribute__((unused)), uint rows){
    return sizeof(Matrix) + rows*(sizeof(Vector*) + VECTOR_CONTAINER_SIZE);
}

Matrix *Matrix_slice(Matrix *M, char *memory, uint row_min, uint row_max, uint column_min, uint column_max){
    assert(row_min <= row_max && row_max <= M->rows);
    assert(column_min <= column_max && column_max <= M->columns);
    Matrix *result = (Matrix*)memory;
    uint num_rows = row_max - row_min;
    uint num_cols = column_max - column_min;
    result->p = M->p;
    result->rows = num_rows;
    result->columns = num_cols;
    result->vectors = (Vector**)(result + 1);
    Vector **matrix_ptr = (Vector**)result->vectors; 
    char *vector_ptr = (char*)(matrix_ptr + num_rows);
    Vector *initialized_vector_ptr;
    for(uint i = 0; i < num_rows; i++){
        initialized_vector_ptr = Vector_initialize(M->p, &vector_ptr, 0, 0);
        *matrix_ptr = initialized_vector_ptr;
        Vector_slice(initialized_vector_ptr, M->vectors[i], column_min, column_max);
        matrix_ptr ++;
    }
    assert(matrix_ptr == (Vector **)(result->vectors + num_rows));
    assert(vector_ptr == (char*)matrix_ptr + num_rows * VECTOR_CONTAINER_SIZE);
    return result;
}

uint Matrix_toString(char *buffer, Matrix *M){
    int len = 0;
    len += sprintf(buffer + len, "    [\n");
    for(uint i = 0; i < M->rows; i++){
        len += sprintf(buffer + len, "        ");
        len += Vector_toString(buffer + len, M->vectors[i]);
        len += sprintf(buffer + len, ",\n");
    }
    len += sprintf(buffer + len, "    ]\n");
    return len;
}

void Matrix_print(Matrix *matrix){
    char buffer[10000];
    Matrix_toString(buffer, matrix);
    printf("%s\n", buffer);
}

void Matrix_printSlice(Matrix *M, uint col_end, uint col_start){
    for(uint i = 0; i < M->rows; i++){
        char buffer[2000];
        uint len = 0;
        char slice_memory[VECTOR_CONTAINER_SIZE];
        char *slice_ptr = slice_memory;
        Vector *slice = Vector_initialize(M->p, &slice_ptr, 0, 0);
        len += sprintf(buffer + len, "    ");
        Vector_slice(slice, M->vectors[i], 0, col_end);
        len += Vector_toString(buffer + len, slice);
        len += sprintf(buffer + len, "; ");
        Vector_slice(slice, M->vectors[i], col_start, M->columns);
        len += Vector_toString(buffer + len, slice);
        printf("%s\n", buffer);
    }
    printf("\n");
}

void Matrix_getRowPermutation(Matrix *M, uint *result){
    Vector *first_vector = (Vector*)(M->vectors + M->rows);
    for(uint i=0; i < M->rows; i++){
        uint j = ((uint64)M->vectors[i] - (uint64)first_vector)/VECTOR_CONTAINER_SIZE; // why is this sizeof(VectorPrivate)??
        result[i] = j;
    }
}

void Matrix_applyRowPermutation(Matrix *M, uint *permutation, uint rows){
    Vector *temp[rows];
    for(uint i=0; i < rows; i++){
        temp[i] = M->vectors[permutation[i]];
    }
    memcpy(M->vectors, temp, rows * sizeof(Vector*));
}


void rowReduce(Matrix *M, int *column_to_pivot_row, uint col_end, uint col_start){
    Vector **matrix = M->vectors;
    uint p = M->p;
    uint columns = M->columns;
    uint rows = M->rows;
    // Fill matrix with -1s = 0xFFFFFFFF
    memset(column_to_pivot_row, 0xFF, columns * sizeof(int));
    if(rows == 0){
        return;
    }
    VectorIterator rowIterators[rows];
    for(uint i = 0; i < rows; i++){
        rowIterators[i] = Vector_getIterator(matrix[i]);
    }
    uint pivot = 0;
    for(uint pivot_column = 0; pivot_column < columns; pivot_column++){
        // Search down column for a nonzero entry.
        uint pivot_row;
        for(pivot_row = pivot; pivot_row < rows; pivot_row++){
            if(rowIterators[pivot_row].value != 0){
                break;
            }
        }
        // No pivot in pivot_column.
        if(pivot_row == rows){
            for(uint i = 0; i < rows; i++){
                rowIterators[i] = Vector_stepIterator(rowIterators[i]);
            }
            continue;
        }
        // Record position of pivot.
        column_to_pivot_row[pivot_column] = pivot;

        if(col_end > 0){
            Matrix_printSlice(M, col_end, col_start);
        }

        // Pivot_row contains a row with a pivot in current column.
        // Swap pivot row up.
        Vector *temp = matrix[pivot];
        matrix[pivot] = matrix[pivot_row];
        matrix[pivot_row] = temp;

        VectorIterator temp_it = rowIterators[pivot];
        rowIterators[pivot] = rowIterators[pivot_row];
        rowIterators[pivot_row] = temp_it;

        if(col_end > 0){
            printf("row(%d) <==> row(%d)\n", pivot, pivot_row);
            Matrix_printSlice(M, col_end, col_start);
        }

        // Divide pivot row by pivot entry
        int c = rowIterators[pivot].value;
        int c_inv = inverse(p, c);
        Vector_scale(matrix[pivot], c_inv);

        if(col_end > 0){
            printf("row(%d) *= %d\n", pivot, c_inv);
            Matrix_printSlice(M, col_end, col_start);
        }
        for(uint i = 0; i < rows; i++){
            // Between pivot and pivot_row, we already checked that the pivot column is 0, so skip ahead a bit.
            if(i == pivot){
                i = pivot_row;
                continue;
            }
            int pivot_column_entry = rowIterators[i].value;
            int row_op_coeff = (p - pivot_column_entry) % p;
            if(row_op_coeff == 0){
                continue;
            }
            // Do row operation
            Vector_add(matrix[i], matrix[pivot], row_op_coeff);
            if(col_end > 0){
                printf("row(%d) <== row(%d) + %d * row(%d)\n", i, i, row_op_coeff, pivot);
                Matrix_printSlice(M, col_end, col_start);
            }
        }
        pivot ++;
        for(uint i = 0; i < rows; i++){
            rowIterators[i] = Vector_stepIterator(rowIterators[i]);
        }        
    }
    return;
}



Subspace *Subspace_construct(uint p, uint rows, uint columns){
    assert(columns < MAX_DIMENSION);
    Subspace *k = malloc(
        sizeof(Subspace) 
        + columns * sizeof(uint)
        + Matrix_getSize(p, rows, columns)
    );
    k->column_to_pivot_row = (int*)(k + 1);
    k->matrix = Matrix_initialize((char*)(k->column_to_pivot_row + columns), p, rows, columns);
    return k;
}

void Subspace_serialize(char **buffer, Subspace *kernel){
    memcpy(*buffer, kernel, sizeof(*kernel));
    *buffer += sizeof(Subspace);
    Matrix_serialize(buffer, kernel->matrix);    
    memcpy(*buffer, kernel->column_to_pivot_row, kernel->matrix->columns * sizeof(int));
    *buffer += kernel->matrix->columns * sizeof(int);
}

Subspace *Subspace_deserialize(char **buffer){
    Subspace *result = (Subspace*) *buffer;
    *buffer += sizeof(Subspace);
    result->matrix = Matrix_deserialize(buffer);
    result->column_to_pivot_row = (int*) *buffer;
    *buffer += result->matrix->columns * sizeof(int);
    return result;
}

size_t Subspace_getSize(uint p, uint rows, uint columns){
    return sizeof(Subspace) + columns * sizeof(int) + Matrix_getSize(p, rows, columns);
}

void Subspace_free(Subspace *k){
    free(k);
}

// Calling code responsible for freeing returned Subspace.
// matrix -- a row reduced augmented matrix
// column_to_pivot_row -- the pivots in matrix (also returned by row_reduce)
// first_source_column -- which block of the matrix is the source of the map
Subspace *Matrix_computeKernel(Matrix *matrix, int *column_to_pivot_row, uint first_source_column){
    const uint p = matrix->p;
    const uint source_dimension = matrix->columns - first_source_column;

    // Find the first kernel row
    uint first_kernel_row = matrix->rows;
    for(uint i = first_source_column; i < matrix->columns; i ++){
        if(column_to_pivot_row[i] != -1){
            first_kernel_row = column_to_pivot_row[i];
            break;
        }
    }
    // Every row after the first kernel row is also a kernel row, so now we know how big it is and can allocate space.
    uint kernel_dimension = matrix->rows - first_kernel_row;
    Subspace *kernel = Subspace_construct(p, kernel_dimension, source_dimension);
    if(kernel_dimension == 0){
        memset(kernel->column_to_pivot_row, 0xFF, kernel->matrix->columns * sizeof(int));
        return kernel;
    }
    // Write pivots into kernel
    for(uint i = 0; i < source_dimension; i++){
        // Turns -1 into some negative number... make sure to check <0 for no pivot in column...
        kernel->column_to_pivot_row[i] = column_to_pivot_row[i + first_source_column] - first_kernel_row;
    }
    // array_print("col_to_piv_row: %s\n", (uint*)kernel->column_to_pivot_row, kernel->matrix->columns);

    // Copy kernel matrix into kernel
    for(uint row = 0; row < kernel_dimension; row++){
        char slice_memory[Vector_getSize(p, 0, 0)];
        char *slice_ptr = slice_memory;
        Vector *slice = Vector_initialize(p, &slice_ptr, 0, 0);
        Vector_slice(slice, matrix->vectors[first_kernel_row + row], first_source_column, first_source_column + source_dimension);
        Vector_assign(kernel->matrix->vectors[row], slice);
    }
    return kernel;
}


uint Matrix_extendImage(Matrix *matrix, uint first_source_column, uint source_dimension, uint first_empty_row, int *current_pivots, Subspace *desired_image){
    uint p = matrix->p;
    uint homology_dimension = 0;
    int *desired_pivots = desired_image->column_to_pivot_row;
    for(uint i = 0; i < desired_image->matrix->columns; i++){
        assert(current_pivots[i] < 0 || desired_pivots[i] >= 0);
        if(current_pivots[i] < 0 && desired_pivots[i] >= 0){
            // Look up the cycle that we're missing and add a generator hitting it.
            int kernel_vector_row = desired_pivots[i];
            // assert(kernel_vector_row < previous_kernel->kernel->rows);
            Vector *new_image = desired_image->matrix->vectors[kernel_vector_row];
            // Write new image to full_matrix
            Vector *matrix_row = matrix->vectors[first_empty_row];
            // Write elementary basis vector into source block of full_matrix
            Vector_setToZero(matrix_row);
            Vector_setEntry(matrix_row, first_source_column + source_dimension + homology_dimension, 1);
            // Stack allocate slice
            char slice_memory[Vector_getSize(p, 0, 0)];
            char *slice_ptr = slice_memory;
            Vector *slice = Vector_initialize(p, &slice_ptr, 0, 0);
            Vector_slice(slice, matrix_row, 0, desired_image->matrix->columns);
            // Write new_image into slice.
            Vector_assign(slice, new_image);
            first_empty_row++;
            homology_dimension++;
        }
    }
    return homology_dimension;
}

void Matrix_quasiInverse_apply(Vector *target, Subspace *image, Matrix *quasi_inverse, Vector *input){
    uint row = 0;
    for(uint i = 0; i < image->matrix->columns; i++){
        if(image->column_to_pivot_row[i] < 0){
            continue;
        }
        uint coeff = Vector_getEntry(input, i);
        Vector_add(target, quasi_inverse->vectors[row], coeff);
        row ++;
    }
}