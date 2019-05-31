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

void initializeModpLookupTable(uint p);
void initializeLimbBitIndexLookupTable(uint p);

uint modPLookup(uint p, uint n);

typedef struct VectorInterface_s VectorInterface;

typedef struct  {
    struct VectorInterface_s *interface;
    uint p;
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

struct VectorInterface_s {    
    uint container_size;
    uint (*getEntriesPer64Bits)(uint p);
    size_t (*getSize)(uint p, uint dimension, uint offset);
    
    Vector *(*initialize)(uint p, char *vector_container, char *memory, uint dimension, uint offset);
    Vector *(*construct)(uint p, uint dimension, uint offset);
    void (*free)(Vector *v); 

    void (*assign)(Vector *target, Vector *source);
    void (*setToZero)(Vector *target);
    void (*pack)(Vector *target, uint *source);
    void (*unpack)(uint *target, Vector *source);

    uint (*getEntry)(Vector *v, uint index);
    void (*setEntry)(Vector *v, uint index, uint value);
    
    void (*slice)(Vector *result, Vector *source, uint start, uint end);
    VectorIterator (*getIterator)(Vector *v); 
    VectorIterator (*stepIterator)(VectorIterator);

    void (*addBasisElement)(Vector *target, uint idx, uint c);
    void (*addArray)(Vector *target, uint *source, uint c);
    void (*add)(Vector *target, Vector *source, uint c);
    void (*scale)(Vector *target, uint c);   
};

extern VectorInterface VectorGenericInterface;
extern VectorInterface Vector2Interface;

uint getEntriesPer64Bits(uint p);
size_t Vector_getSize(uint p, uint dimension, uint offset);

Vector *Vector_initialize(VectorInterface *interface, uint p, char *vector_container, char *memory, uint dimension, uint offset);
Vector *Vector_construct(VectorInterface *interface, uint p, uint dimension, uint offset);

Vector *VectorGeneric_initialize(uint p, char *vector_container, char *memory, uint dimension, uint offset);
Vector *VectorGeneric_construct(uint p, uint dimension, uint offset);

Vector *Vector2_initialize(uint p, char *vector_container, char *memory, uint dimension, uint offset);
Vector *Vector2_construct(uint p, uint dimension, uint offset);

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

// Maybe delete addArray methods?
void VectorGeneric_addBasisElement(Vector *target, uint idx, uint c);
void VectorGeneric_addArray(Vector *target, uint *source, uint c);
void VectorGeneric_add(Vector *target, Vector *source, uint c);
void VectorGeneric_scale(Vector *target, uint c);

void Vector2_addBasisElement(Vector *target, uint idx, uint c);
void Vector2_addArray(Vector *target, uint *source, uint c);
void Vector2_add(Vector *target, Vector *source, uint c);
void Vector2_scale(Vector *target, uint c);


uint Vector_toString(char *buffer, Vector *v);
void Vector_print(Vector *v);

typedef struct {
    uint p;
    uint rows;
    uint columns;
    Vector **matrix;
} Matrix;

uint Matrix_getSize(VectorInterface *interface, uint p, uint rows, uint cols);
Matrix *Matrix_initialize(char *memory, VectorInterface *interface, uint p, uint rows, uint cols);
Matrix *Matrix_construct(VectorInterface *interface, uint p, uint rows, uint cols);

Matrix *MatrixGeneric_construct(uint p, uint rows, uint cols);
Matrix *Matrix2_construct(uint p, uint rows, uint cols);

uint Matrix_toString(char *buffer, Matrix *M);
void Matrix_print(Matrix *matrix);
uint Matrix_getSliceSize(uint rows);
Matrix *Matrix_slice(Matrix *M, char *memory, uint row_min, uint row_max, uint column_min, uint column_max);
void Matrix_printSlice(Matrix *M, uint col_end, uint col_start);

void rowReduce(Matrix *M, int *column_to_pivot_row, uint, uint);

#endif //C_FPVECTOR_H
