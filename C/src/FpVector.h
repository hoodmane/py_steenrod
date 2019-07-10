//
// Created by Hood on 5/22/2019.
//

/**
 * This file defines our vector container type. The interface defines basic operations
 * like add, scale, assign, slice, iterate. 
 * 
 * This interface leaves us the option to do something 
 * more efficient later if possible (graphics card vector operations? different bit 
 * packing scheme?) 
 * 
 * The current representation: when p = 2, each bit is an entry and we use bit xor to add.
 * When p is odd, bit pack them so that we have up to p*(p-1) values we can store without overflow.
 * Then we can do v <- v + c*w all at once without overflow and split it into chunks and reduce when we're done.
 * I'm not sure this is better than packing it as densely as possible -- it'd be a reasonable experiment
 * to compare the options.
 * There seemed to be a 20% performance gain to this bitpacking scheme at the prime 3
 * as compared to the trivial strategy of using a char* when I tested it. The perfomance 
 * improvement is surely even larger at the prime 2, but that's a bit unfair.
 * 
 * We use the pointer cast to have private fields. When stack allocating vectors, calling code
 * is expected to ask FpVector how big it is by calling Vector_getContainerSize(p).
 */


#ifndef C_FP_VECTOR_H
#define C_FP_VECTOR_H

#include <stdbool.h>
#include <stdlib.h>
// The code will always segfault if this is 150000, so this is about as much memory
// as I can ever stack allocate.
// TODO: Move this stuff to global variables.
#define MAX_DIMENSION 147500
typedef unsigned long long uint64;
typedef long long int64;
typedef unsigned int uint;

// Write a string representation of the array A of length length into buffer.
int array_toString(char *buffer, uint *A, uint length);

// Print array A of length length into %s in format string (should only contain a single %s).
void array_print(char *format_string, uint *A, uint length);

// Look up n % p in a table. Not sure if this helps at all. If it does help it should be moved to combinatorics
uint modPLookup(uint p, uint n); 

// This is the public Vector interface. The backing is hidden after this struct.
typedef struct  {
    uint dimension; // Dimension of space
    uint size;   // How big is it actually?
    uint offset; // offset is zero unless we took a slice of something. For handling slices.
                 // offset is a number of bits.
} Vector;

typedef struct {
    bool has_more;
    uint index;
    uint value;
// private fields 
    Vector *vector;
    uint limb_index;
    uint bit_index;    
} VectorIterator;


// Get the size of the backing of a Vector. 
// Used for stack allocating vectors or writing them into tables 
// (Vector_construct does a separate call to malloc for each vector allocated.)
size_t Vector_getSize(uint p, uint dimension, uint offset);

// This gets rounds the dimension up to be divisible by EntriesPer64Bits
// so that if we slice the next chunk of the vector, it will have offset 0.
// Used for making block vectors.
uint Vector_getPaddedDimension(uint p, uint dimension, uint offset);

/**
 * Since we use a lot of short lived vectors, we want to be able to stack allocate them.
 * We also want to allocate large tables.
 * To stack allocate a vector:
 *      size_t container_size = Vector_getContainerSize(p);
 *      size_t total_size = container_size + Vector_getSize(p, dimension, offset);
 *      char memory[total_size];
 *      myVector = Vector_initialize(p, memory, memory + container_size, dimension, offset);
 */
Vector *Vector_initialize(uint p, char *memory, uint dimension, uint offset);
Vector *Vector_construct(uint p, uint dimension, uint offset);

void Vector_free(Vector *v); 

// Vector_assign and Vector_setToZero are unsuitable for use on a slice
// unless that slice is a "block" -- that is, it should start in a position with offset = 0
// and either go to the end of the source vector or the source vector should be padded
// with garbage fields up until the end of the next leg.
void Vector_assign(Vector *target, Vector *source);
void Vector_setToZero(Vector *target);

// Convert between arrays and vectors.
void Vector_pack(Vector *target, uint *source);
void Vector_unpack(uint *target, Vector *source);

uint Vector_getEntry(Vector *v, uint index);
void Vector_setEntry(Vector *v, uint index, uint value);


/**
 * To slice a vector:
 *      char slice_memory[Vector_getContainerSize(p)];
 *      Vector *slice = Vector_initialize(p, slice_memory, NULL, 0, 0);
 *      Vector_slice(slice, source, start, end);
 * The slice occupies the same memory as the original vector -- this is NOT A COPY!
 */
void Vector_slice(Vector *result, Vector *source, uint start, uint end);

/**
 * To iterate over a vector:
 *      for(
 *          VectorIterator it = Vector_getIterator(myVector);
 *          it.has_more;
 *          it = Vector_stepIterator(it)
 *      ){
 *          // Current entry is it.value, current index is it.index.
 *      }
 * I don't think this provides any performance benefit over using GetEntry, but let's pretend that it does.
 */
VectorIterator Vector_getIterator(Vector *v); 
VectorIterator Vector_stepIterator(VectorIterator);

// Add c to entry at idx.
void Vector_addBasisElement(Vector *target, uint idx, uint c);
// Add c*source[i] to entry i of target. source should be at least as long as target->dimension
// This is a kind of dumb function, but it's currently used in MilnorAlgebra_multiply. We don't bother optimizing it.
void Vector_addArray(Vector *target, uint *source, uint c);
// Add target <- target + c*source. target and source should be the same dimension and offset.
// If source is bigger than target, you need to slice out a chunk of the source first and vice versa.
// Offsets take some care to deal with.
void Vector_add(Vector *target, Vector *source, uint c);
void Vector_scale(Vector *target, uint c);

// Convert vector to a string and write into buffer
uint Vector_toString(char *buffer, Vector *v);
void Vector_print(char *fmt_string, Vector *v);

uint Vector_serialize(char *buffer, Vector *v);

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

uint Matrix_toString(char *buffer, Matrix *M);
void Matrix_print(Matrix *matrix);
uint Matrix_getSliceSize(uint p, uint rows);
Matrix *Matrix_slice(Matrix *M, char *memory, uint row_min, uint row_max, uint column_min, uint column_max);
void Matrix_printSlice(Matrix *M, uint col_end, uint col_start);

// Row reduce M. For each column i, column_to_pivot_row[i] is equal to: the row
// with a pivot in column i or -1 if no such row exists.
// When you delete the -1's from column_to_pivot_row, it looks like [0,1,2,...,rank(M)].
void rowReduce(Matrix *M, int *column_to_pivot_row, uint, uint);


typedef struct {
    int *column_to_pivot_row;
    Matrix *kernel;
} Kernel;

Kernel *Kernel_construct(uint p, uint rows, uint columns);
void Kernel_free(Kernel *k);

#endif //C_FP_VECTOR_H
