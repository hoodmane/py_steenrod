//
// Created by Hood on 5/22/2019.
//

#ifndef C_FPVECTOR_H
#define C_FPVECTOR_H

#include <stdbool.h>
#include <stdlib.h>
#define MAX_DIMENSION 100000
typedef unsigned long long uint64;
typedef long long int64;
typedef unsigned int uint;
int array_toString(char *buffer, uint *A, uint length);
void array_print(uint *A, uint length);


uint modPLookup(uint p, uint n);

typedef struct  {
    uint dimension;
    uint size;
    uint offset;
} Vector;

typedef struct {
    bool has_more;
    uint index;
    uint value;
// private
    Vector *vector;
    uint limb_index;
    uint bit_index;    
} VectorIterator;

uint Vector_getContainerSize(uint p);
size_t Vector_getSize(uint p, uint dimension, uint offset);
uint Vector_getPaddedDimension(uint p, uint dimension, uint offset);

Vector *Vector_initialize(uint p, char *vector_container, char *memory, uint dimension, uint offset);
Vector *Vector_construct(uint p, uint dimension, uint offset);

void Vector_free(Vector *v); 

void Vector_assign(Vector *target, Vector *source);
void Vector_setToZero(Vector *target);
void Vector_pack(Vector *target, uint *source);
void Vector_unpack(uint *target, Vector *source);

uint Vector_getEntry(Vector *v, uint index);
void Vector_setEntry(Vector *v, uint index, uint value);

void Vector_slice(Vector *result, Vector *source, uint start, uint end);
VectorIterator Vector_getIterator(Vector *v); 
VectorIterator Vector_stepIterator(VectorIterator);

void Vector_addBasisElement(Vector *target, uint idx, uint c);
void Vector_addArray(Vector *target, uint *source, uint c);
void Vector_add(Vector *target, Vector *source, uint c);
void Vector_scale(Vector *target, uint c);

uint Vector_toString(char *buffer, Vector *v);
void Vector_print(Vector *v);

typedef struct {
    uint p;
    uint rows;
    uint columns;
    Vector **matrix;
} Matrix;

uint Matrix_getSize(uint p, uint rows, uint cols);
Matrix *Matrix_initialize(char *memory, uint p, uint rows, uint cols);
Matrix *Matrix_construct(uint p, uint rows, uint cols);
void Matrix_free(Matrix *M);

Matrix *MatrixGeneric_construct(uint p, uint rows, uint cols);
Matrix *Matrix2_construct(uint p, uint rows, uint cols);

uint Matrix_toString(char *buffer, Matrix *M);
void Matrix_print(Matrix *matrix);
uint Matrix_getSliceSize(uint p, uint rows);
Matrix *Matrix_slice(Matrix *M, char *memory, uint row_min, uint row_max, uint column_min, uint column_max);
void Matrix_printSlice(Matrix *M, uint col_end, uint col_start);

void rowReduce(Matrix *M, int *column_to_pivot_row, uint, uint);

#endif //C_FPVECTOR_H
