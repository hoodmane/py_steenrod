//
// Created by Hood on 5/22/2019.
//

#include "FpVector.h"
#include "combinatorics.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int array_toString(char *buffer, uint *A, uint length){
    buffer[0] = '[';
    buffer[1] = '\0';
    int len = 1;
    for (int i = 0; i < length; i++) {
        len += sprintf(buffer + len, "%d, ", A[i]);
    }
    len += sprintf(buffer + len, "]");
    return len;
}

void array_print(uint *A, uint length){
    char buffer[1000];
    array_toString(buffer, A, length);
    printf("%s", buffer);
}


typedef struct VectorImplementation_s VectorImplementation;
VectorImplementation *getVectorImplementation(uint p);

typedef struct {
    uint dimension;
    uint size;
    uint offset;
// Private fields:
    struct VectorImplementation_s *implementation;
    uint number_of_limbs;
    uint64 *vector;
} VectorStd;

struct VectorImplementation_s {  
    uint p;  
    void (*addBasisElement)(Vector *target, uint idx, uint c);
    void (*addArray)(Vector *target, uint *source, uint c);
    void (*add)(Vector *target, Vector *source, uint c);
    void (*scale)(Vector *target, uint c);   
};

// Generated with Mathematica:
//     Ceiling[Log2[# (# - 1) + 1 &[Prime[Range[54]]]]]
// But for 2 it should be 1.
uint bitlengths[MAX_PRIME_INDEX] = { 
     1, 3, 5, 6, 7, 8, 9, 9, 9, 10, 10, 11, 11, 11, 12, 12, 12, 12, 13,     
     13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15,    
     15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16 
};

uint getBitlength(uint p){
    return bitlengths[prime_to_index_map[p]];
}

// Generated with Mathematica:
//     2^Ceiling[Log2[# (# - 1) + 1 &[Prime[Range[54]]]]]-1
// But for 2 it should be 1.
uint bitmasks[MAX_PRIME_INDEX] = {
    1, 7, 31, 63, 127, 255, 511, 511, 511, 1023, 1023, 2047, 2047, 2047, 
    4095, 4095, 4095, 4095, 8191, 8191, 8191, 8191, 8191, 8191, 16383, 
    16383, 16383, 16383, 16383, 16383, 16383, 32767, 32767, 32767, 32767, 
    32767, 32767, 32767, 32767, 32767, 32767, 32767, 65535, 65535, 65535,
    65535, 65535, 65535, 65535, 65535, 65535, 65535, 65535, 65535
};

uint getBitMask(uint p){
    return bitmasks[prime_to_index_map[p]];
}

// Generated with Mathematica:
//      Floor[64/Ceiling[Log2[# (# - 1) + 1 &[Prime[Range[54]]]]]]
// But for 2 it should be 64.
uint entries_per_64_bits[MAX_PRIME_INDEX] = {
    64, 21, 12, 10, 9, 8, 7, 7, 7, 6, 6, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4,  
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

uint getEntriesPer64Bits(uint p){
    return entries_per_64_bits[prime_to_index_map[p]];   
}

uint *modplookuptable[MAX_PRIME_INDEX] = {0};

// Called by initializePrime
void initializeModpLookupTable(uint p){
    uint p_times_p_minus_1 = p*(p-1) + 1;
    uint *table = malloc((p_times_p_minus_1 + 1) * sizeof(uint));
    for(uint i = 0; i <= p_times_p_minus_1; i++){
        table[i] = i % p;
    }
    modplookuptable[prime_to_index_map[p]] = table;
}

// n must be in the range 0 <= n <= p * (p-1)
uint modPLookup(uint p, uint n){
    // printf("%d ", n);
    return modplookuptable[prime_to_index_map[p]][n];
}

uint Vector_getContainerSize(uint p){
    return sizeof(VectorStd);
}

typedef struct {
    uint limb;
    uint bit_index;
} LimbBitIndexPair;

LimbBitIndexPair *limbBitIndexLookupTable[MAX_PRIME_INDEX] = {0};

// Called by initializePrime
void initializeLimbBitIndexLookupTable(uint p){
    uint p_idx = prime_to_index_map[p];
    uint entries_per_limb = getEntriesPer64Bits(p);
    uint bit_length = getBitlength(p);
    LimbBitIndexPair *table = malloc(MAX_DIMENSION * sizeof(LimbBitIndexPair));
    for(uint i = 0; i < MAX_DIMENSION; i++){
        table[i].limb = i/entries_per_limb;
        table[i].bit_index = (i % entries_per_limb) * bit_length;
    }
    limbBitIndexLookupTable[p_idx] = table;
}

LimbBitIndexPair getLimbBitIndexPair(uint p, uint idx){
    assert(prime_to_index_map[p] != -1);
    assert(idx < MAX_DIMENSION);
    return limbBitIndexLookupTable[prime_to_index_map[p]][idx];
}

size_t Vector_getSize(uint p, uint dimension, uint offset){
    assert(dimension < MAX_DIMENSION);
    assert(offset < 64);
    uint bit_length = getBitlength(p);
    size_t size = (dimension == 0) ? 0 : (getLimbBitIndexPair(p, dimension + offset/bit_length - 1).limb + 1);
    size *= sizeof(uint64);
    return size; // Need extra space for vector fields.
}

uint Vector_getPaddedDimension(uint p, uint dimension, uint offset){
    uint entries_per_limb = getEntriesPer64Bits(p);
    return ((dimension + offset + entries_per_limb - 1)/entries_per_limb)*entries_per_limb;
}


Vector *Vector_initialize(uint p, char *container, char *memory, uint dimension, uint offset){
    VectorImplementation *vectImpl = getVectorImplementation(p);
    VectorStd *v = (VectorStd *) container;
    v->implementation = vectImpl;
    v->dimension = dimension;
    v->size = Vector_getSize(p, dimension, offset);
    v->offset = offset;
    v->number_of_limbs = (v->size + sizeof(uint64) - 1)/sizeof(uint64);
    v->vector = dimension == 0 ? NULL : (uint64*)memory;
    memset(v->vector, 0, v->size);
    return (Vector*)v;
}

// There is no case distinction between Vector2 and VectorGeneric for the functions that just 
// get and set values.
Vector *Vector_construct(uint p, uint dimension, uint offset){
    uint size = Vector_getSize(p, dimension, offset);
    size += Vector_getContainerSize(p);
    char *memory = malloc(size);
    Vector *result = Vector_initialize(p, memory, (char*)((VectorStd*)memory + 1), dimension, offset);
    return result;
}

void Vector_free(Vector *vector){
    free(vector);
}

void Vector_setToZero(Vector *target){
    assert(target->offset == 0); // setToZero doesn't handle slices right now.
    VectorStd *t = (VectorStd*) target;
    memset(t->vector, 0, t->size);
}


void Vector_assign(Vector *target, Vector *source){
    assert(source->dimension == target->dimension);
    assert(source->offset == target->offset);
    assert(source->offset == 0); // Vector_assign doesn't handle slices right now.
    VectorStd *t = (VectorStd*) target;
    VectorStd *s = (VectorStd*) source;
    memcpy(t->vector, s->vector, s->size);
}

uint unpackLimbHelper(uint *limb_array, VectorStd *vector, uint limb_idx, uint bit_min, uint bit_max){
    uint bit_mask = getBitMask(vector->implementation->p);
    uint bit_length = getBitlength(vector->implementation->p);
    uint64 limb_value = vector->vector[limb_idx];
    uint idx = 0;
    for(uint j = bit_min; j < bit_max - bit_length + 1; j += bit_length){
        limb_array[idx] = (limb_value >> j) & bit_mask;
        idx++;
    }
    return idx;
}

uint unpackLimb(uint *limb_array, VectorStd *vector, uint limb_idx){
    uint bit_length = getBitlength(vector->implementation->p);
    uint entries_per_64_bits = getEntriesPer64Bits(vector->implementation->p);
    uint bit_min = 0;
    uint bit_max = 64;    
    if(limb_idx == 0){
        bit_min = vector->offset;
    }
    if(limb_idx == vector->number_of_limbs - 1){
        bit_max = (vector->offset + vector->dimension * bit_length)%(bit_length * entries_per_64_bits);;
        if(bit_max == 0){
            bit_max = 64;
        }    
    }
    return unpackLimbHelper(limb_array, vector, limb_idx, bit_min, bit_max);
}

uint packLimbHelper(VectorStd *vector, uint *limb_array, uint limb_idx, uint bit_min, uint bit_max, uint64 bit_mask){
    uint bit_length = getBitlength(vector->implementation->p);
    uint idx = 0;
    uint64 limb_value = vector->vector[limb_idx] & bit_mask;
    for(uint j = bit_min; j < bit_max - bit_length + 1; j += bit_length){
        limb_value |= ((uint64) limb_array[idx]) << j;
        idx ++;
    }
    vector->vector[limb_idx] = limb_value;
    return idx;
}

uint packLimb(VectorStd *vector, uint *limb_array, uint limb_idx){
    uint bit_length = getBitlength(vector->implementation->p);
    uint entries_per_64_bits = getEntriesPer64Bits(vector->implementation->p);
    uint bit_min = 0;
    uint bit_max = 64;    
    if(limb_idx == 0){
        bit_min = vector->offset;
    }
    if(limb_idx == vector->number_of_limbs - 1){
        bit_max = (vector->offset + vector->dimension * bit_length)%(bit_length * entries_per_64_bits);
        if(bit_max == 0){
            bit_max = 64;
        }    
    }
    uint64 bit_mask = 0;
    if(bit_max - bit_min < 64){
        bit_mask = (1LL << (bit_max - bit_min)) - 1;
        bit_mask <<= bit_min;
        bit_mask = ~bit_mask;
    }
    return packLimbHelper(vector, limb_array, limb_idx, bit_min, bit_max, bit_mask);
}


void Vector_pack(Vector *target, uint *source){
    assert(target != NULL);
    VectorStd * t = (VectorStd*) target;
    for(uint limb = 0; limb < t->number_of_limbs; limb++){
        source += packLimb(t, source, limb);
    }
}

void Vector_unpack(uint * target, Vector * source){ 
    assert(target != NULL);
    VectorStd * s = (VectorStd*) source;    
    for(uint i = 0; i < s->number_of_limbs; i++){
        target += unpackLimb(target, s, i);
    }
}

uint Vector_getEntry(Vector *vector, uint index){
    assert(vector != NULL);
    assert(index < vector->dimension);
    VectorStd *v = (VectorStd*) vector;    
    uint64 bit_mask = getBitMask(v->implementation->p);
    LimbBitIndexPair limb_index = getLimbBitIndexPair(v->implementation->p, index + v->offset);
    uint64 result = v->vector[limb_index.limb];
    result >>= limb_index.bit_index;
    result &= bit_mask;
    return result;
}

void Vector_setEntry(Vector *vector, uint index, uint value){
    assert(vector != NULL);
    assert(index < vector->dimension);
    VectorStd *v = (VectorStd*) vector;    
    uint64 bit_mask = getBitMask(v->implementation->p);
    LimbBitIndexPair limb_index = getLimbBitIndexPair(v->implementation->p, index + v->offset);
    uint64 *result = &(v->vector[limb_index.limb]);
    *result &= ~(bit_mask << limb_index.bit_index);
    *result |= (((uint64)value) << limb_index.bit_index);
}

void Vector_slice(Vector *result, Vector *source, uint start, uint end){
    assert(start <= end);
    assert(end <= source->dimension);
    VectorStd *r = (VectorStd*) result;
    VectorStd *s = (VectorStd*) source;
    r->implementation = s->implementation;
    r->dimension = end - start;
    if(start == end){
        r->size = 0;
        r->number_of_limbs = 0;
        r->vector = NULL;
        return;    
    }
    LimbBitIndexPair limb_index = getLimbBitIndexPair(r->implementation->p, start + source->offset);
    r->offset = limb_index.bit_index;
    r->size = Vector_getSize(r->implementation->p, r->dimension, r->offset);
    r->number_of_limbs = r->size/sizeof(uint64);
    r->vector = s->vector + limb_index.limb;
}

VectorIterator Vector_getIterator(Vector *vector){
    VectorIterator result;
    VectorStd *v = (VectorStd*)vector;
    if(v->dimension == 0){
        result.has_more = false;
        return result;
    }
    result.has_more = true;
    result.vector = vector;
    result.index = 0;
    result.limb_index = 0;
    result.bit_index = v->offset;
    uint bit_mask = getBitMask(v->implementation->p);
    result.value = ((v->vector[result.limb_index]) >> result.bit_index) & bit_mask;
    return result;
}

VectorIterator Vector_stepIterator(VectorIterator it){
    it.index ++;
    VectorStd *v = (VectorStd*)it.vector;
    it.has_more = it.index < it.vector->dimension;
    if(!it.has_more){
        return it;
    }
    uint bit_mask = getBitMask(v->implementation->p);
    uint bit_length = getBitlength(v->implementation->p);
    it.bit_index += bit_length;
    if(it.bit_index >= 64 - bit_length + 1){
        it.limb_index ++;
        it.bit_index = 0;
    }
    it.value = ((v->vector[it.limb_index]) >> it.bit_index) & bit_mask;
    return it;
}

// For the arithmetic, at 2 we xor whereas at odd primes we have to do a bit more work.
// Odd prime vector arithmetic
// We've chosen our packing so that we have enough space to fit p*(p-1) in each limb.
// So we do the arithmetic in place. Then we have to reduce each slot by p.
// We do this by pulling the slots out and reducing each one separately mod p via table lookup.

void VectorGeneric_addBasisElement(Vector *target, uint index, uint coeff){
    assert(index < target->dimension);
    VectorStd *t = (VectorStd*) target;    
    uint64 bit_mask = getBitMask(t->implementation->p);
    LimbBitIndexPair limb_index = getLimbBitIndexPair(t->implementation->p, index + target->offset);
    uint64 *result = &(t->vector[limb_index.limb]);
    uint64 new_entry = *result >> limb_index.bit_index;
    new_entry &= bit_mask;
    new_entry += coeff;
    new_entry = modPLookup(t->implementation->p, new_entry);
    *result &= ~(bit_mask << limb_index.bit_index);
    *result |= (new_entry << limb_index.bit_index);
}

void VectorGeneric_addArray(Vector *target, uint *source, uint c){
    VectorStd *t = (VectorStd*) target;
    uint source_idx = 0;
    uint entries[getEntriesPer64Bits(t->implementation->p)];   
    for(uint i = 0; i < t->number_of_limbs; i++){
        uint limb_length = unpackLimb(entries, t, i);
        for(uint j = 0; j < limb_length; j++){
            entries[j] = modPLookup(t->implementation->p, entries[j] + c*source[source_idx]);
            source_idx++;
        }
        packLimb(t, entries, i);
    }
}

void VectorGeneric_add(Vector *target, Vector *source, uint coeff){
    assert(source->dimension == target->dimension);
    assert(target->offset == source->offset);
    VectorStd *t = (VectorStd*) target;
    VectorStd *s = (VectorStd*) source;    
    uint entries[getEntriesPer64Bits(t->implementation->p)];    
    for(uint i = 0; i < s->number_of_limbs; i++){
        t->vector[i] = t->vector[i] + coeff * s->vector[i];
        uint limb_length = unpackLimb(entries, t, i);
        for(uint j = 0; j < limb_length; j++){
            entries[j] = modPLookup(t->implementation->p, entries[j]);
        }
        packLimb(t, entries, i);
    }
}

void VectorGeneric_scale(Vector *target, uint coeff){
    // assert(coeff != 0);
    VectorStd *t = (VectorStd*) target;   
    uint entries_per_64_bits = getEntriesPer64Bits(t->implementation->p);    
    uint entries[entries_per_64_bits];
    for(uint i = 0; i < t->number_of_limbs; i++){
        t->vector[i] = coeff * t->vector[i];
        uint limb_length = unpackLimb(entries, t, i);
        for(uint j = 0; j < limb_length; j++){
            entries[j] = modPLookup(t->implementation->p, entries[j]);
        }
        packLimb(t, entries, i);        
    }
}

// Vector2
void Vector2_addBasisElement(Vector *target, uint index, uint coeff){
    assert(index < target->dimension);
    VectorStd *t = (VectorStd*) target;    
    LimbBitIndexPair limb_index = getLimbBitIndexPair(t->implementation->p, index + target->offset);
    uint64 *result = &(t->vector[limb_index.limb]);
    *result ^= ((uint64)coeff << limb_index.bit_index);
}

void Vector2_addArray(Vector *target, uint *source, uint c){
    VectorStd *t = (VectorStd*) target;
    uint bit_length = 1;
    uint source_idx = 0;
    for(uint limb = 0; limb < t->number_of_limbs; limb++){
        uint64 source_limb = 0;
        for(
            uint j = (limb == 0) ? (target->offset) : 0;
            j < (64 - bit_length + 1) && source_idx < target->dimension; 
            j += bit_length
        ){
            source_limb |= ((uint64)source[source_idx] << j);
            source_idx ++;
        }  
        t->vector[limb] ^= source_limb;
    }
}

void Vector2_add(Vector *target, Vector *source, uint coeff){
    assert(
        target->dimension == source->dimension 
        && target->offset == source->offset 
    );    
    VectorStd *t = (VectorStd*) target;
    VectorStd *s = (VectorStd*) source;    
    uint64 *target_ptr = t->vector;
    uint64 *source_ptr = s->vector;
    if(t->number_of_limbs == 0){
        return;
    }
    uint64 source_limb;
    if(t->number_of_limbs == 1){
        source_limb = *source_ptr;
        uint64 bit_mask = -1;
        uint bit_min = source->offset;
        uint bit_max = (source->offset + source->dimension) % 64;
        bit_max = bit_max ? bit_max : 64;
        if(bit_max - bit_min < 64){
            bit_mask = (1LL << (bit_max - bit_min)) - 1;
            bit_mask <<= bit_min;
        }
        source_limb &= bit_mask;
        *target_ptr ^= (coeff * source_limb);
        return;
    }
    source_limb = *source_ptr;    
    // Mask out low order bits
    source_limb &= ~((1LL<<source->offset) - 1);
    *target_ptr ^= coeff * source_limb;    
    target_ptr++;
    source_ptr++;
    for(uint i = 1; i < t->number_of_limbs-1; i++){
        *target_ptr ^= (coeff * (*source_ptr));
        target_ptr++;
        source_ptr++;
    }
    source_limb = *source_ptr;
    // Mask out high order bits
    uint64 bit_mask = -1;
    uint bit_max = (source->offset + source->dimension) % 64;
    if(bit_max){
        bit_mask = (1LL<<bit_max) - 1;
    }
    source_limb &= bit_mask;
    *target_ptr ^= source_limb;
}

// This is a pretty pointless method.
void Vector2_scale(Vector *target, uint coeff){
    assert(coeff != 0);
    return;
    VectorStd *t = (VectorStd*) target;   
    uint64 *target_ptr = t->vector;
    for(uint i = 0; i < t->number_of_limbs; i++){
        *target_ptr *= coeff;
        target_ptr++;
    }
}

VectorImplementation VectorGenericImplementation = {
    0,    
    VectorGeneric_addBasisElement, VectorGeneric_addArray, VectorGeneric_add, VectorGeneric_scale,
};

VectorImplementation Vector2Implementation = {
    0,
    Vector2_addBasisElement, Vector2_addArray, Vector2_add, Vector2_scale,
};

// The generic methods depend on the implementation, so we look it up on the target.
void Vector_addBasisElement(Vector *target, uint idx, uint c){
    ((VectorStd*)target)->implementation->addBasisElement(target, idx, c);
}

void Vector_addArray(Vector *target, uint *source, uint c){
    ((VectorStd*)target)->implementation->addArray(target, source, c);
}
void Vector_add(Vector *target, Vector *source, uint c){
    ((VectorStd*)target)->implementation->add(target, source, c);
}

void Vector_scale(Vector *target, uint c){
    ((VectorStd*)target)->implementation->scale(target, c);
}


VectorImplementation vectorImplementationTable[MAX_PRIME_INDEX];
// Called by initializePrime
void initializeVectorImplementation(uint p){
    if(p == 2){
        vectorImplementationTable[prime_to_index_map[p]] = Vector2Implementation;
    } else {
        vectorImplementationTable[prime_to_index_map[p]] = VectorGenericImplementation;
    }
    vectorImplementationTable[prime_to_index_map[p]].p = p;
}


VectorImplementation *getVectorImplementation(uint p){
    return &vectorImplementationTable[prime_to_index_map[p]];
}

uint Vector_toString(char *buffer, Vector *vector){
    uint len = 0;
    len += sprintf(buffer + len, "[");
    for(VectorIterator it = Vector_getIterator(vector); it.has_more; it = Vector_stepIterator(it)){
        len += sprintf(buffer + len, "%d, ", it.value);
    }
    len += sprintf(buffer + len, "]");
    return len;
}

void Vector_print(Vector *v){
    char buffer[10000];
    Vector_toString(buffer, v);
    printf("%s\n", buffer);
}



uint Matrix_getSize(uint p, uint rows, uint cols){
    // printf("cols: %ud\n", cols);
    assert(cols < MAX_DIMENSION);
    return sizeof(Matrix) 
      + rows * (sizeof(Vector*) + Vector_getContainerSize(p)  + Vector_getSize(p, cols, 0));
}

Matrix *Matrix_initialize(char *memory, uint p, uint rows, uint columns)  {
    uint container_size = Vector_getContainerSize(p);
    uint vector_size = Vector_getSize(p, columns, 0);
    Matrix *matrix = (Matrix*)memory;
    Vector **vector_ptr = (Vector**)(matrix+1);
    char *container_ptr = (char*)(vector_ptr + rows);
    char *values_ptr = container_ptr + rows * container_size;
    matrix->p = p;
    matrix->rows = rows;
    matrix->columns = columns;
    matrix->matrix = vector_ptr;
    for(int row = 0; row < rows; row++){
        *vector_ptr = Vector_initialize(p, container_ptr, values_ptr, columns, 0);
        vector_ptr ++;
        container_ptr += container_size;
        values_ptr += vector_size;
    }
    return matrix;
}

Matrix *Matrix_construct(uint p,  uint rows, uint columns)  {
    char *M = malloc(Matrix_getSize(p, rows, columns));
    // printf("columns: %d, rows: %d\n", columns, rows);
    return Matrix_initialize(M, p, rows, columns);
}

uint Matrix_getSliceSize(uint p, uint rows){
    return sizeof(Matrix) + rows*(sizeof(Vector*) + Vector_getContainerSize(p));
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
    result->matrix = (Vector**)(result + 1);
    VectorStd **matrix_ptr = (VectorStd**)result->matrix; 
    VectorStd *vector_ptr = (VectorStd*)(matrix_ptr + num_rows);
    for(int i = 0; i < num_rows; i++){
        *matrix_ptr = vector_ptr;
        Vector_slice((Vector*)vector_ptr, M->matrix[i], column_min, column_max);
        matrix_ptr ++;
        vector_ptr++;
    }
    assert(matrix_ptr == (VectorStd **)(result->matrix + num_rows));
    assert(vector_ptr == (VectorStd*)matrix_ptr + num_rows);
    return result;
}

uint Matrix_toString(char *buffer, Matrix *M){
    int len = 0;
    len += sprintf(buffer + len, "    [\n");
    for(int i = 0; i < M->rows; i++){
        len += sprintf(buffer + len, "        ");
        len += Vector_toString(buffer + len, M->matrix[i]);
        len += sprintf(buffer + len, ",\n");
    }
    len += sprintf(buffer + len, "    ]\n");
    return len;
}

void Matrix_printMatrix(Matrix *matrix){
    char buffer[10000];
    Matrix_toString(buffer, matrix);
    printf("%s\n", buffer);
}

void Matrix_printSlice(Matrix *M, uint col_end, uint col_start){
    for(uint i = 0; i < M->rows; i++){
        char buffer[2000];
        uint len = 0;
        char slice_memory[Vector_getContainerSize(M->p)];
        Vector *slice = (Vector*)slice_memory;
        len += sprintf(buffer + len, "    ");
        Vector_slice(slice, M->matrix[i], 0, col_end);
        len += Vector_toString(buffer + len, slice);
        len += sprintf(buffer + len, "; ");
        Vector_slice(slice, M->matrix[i], col_start, M->columns);
        len += Vector_toString(buffer + len, slice);
        printf("%s\n", buffer);
    }
    printf("\n");
}

void rowReduce(Matrix *M, int *column_to_pivot_row, uint col_end, uint col_start){
    Vector **matrix = M->matrix;
    uint p = M->p;
    uint columns = M->columns;
    uint rows = M->rows;
    memset(column_to_pivot_row, -1, columns * sizeof(int));
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
        for(int i = 0; i < rows; i++){
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


/**
int main(int argc, char *argv[]){
    // uint p = 3;
    // uint bitmask = 3;
    // uint dim = 49;
    // uint num_vectors = 100;
    // uint scale = 3;    
    // uint num_repeats = 400;
    
    // initializePrime(p);
    // uint arrays[num_vectors][dim];
    // uint container_size = VectorGenericInterface.container_size;
    // uint contents_size = VectorGenericInterface.getSize(p, dim, 0);
    // uint64 vector_container_memory[container_size * num_vectors];
    // uint64 vector_contents_memory[contents_size * num_vectors];
    // uint64 *vector_container_ptr = vector_container_memory;
    // uint64 *vector_contents_ptr = vector_contents_memory;
    // Vector *vectors[num_vectors];
    // for(uint i = 0; i < num_vectors; i++){
    //     for(uint j = 0; j < dim; j++){
    //         arrays[i][j] = modPLookup(p, rand() & bitmask);
    //     }
    // }
    // for(uint i = 0; i < num_vectors; i++){
    //     vectors[i] = VectorGenericInterface.initialize(p, vector_container_ptr, vector_contents_ptr, dim, 0);
    //     vector_container_ptr += container_size;
    //     vector_contents_ptr += contents_size;
    // }
    
    // // for(uint i = 0; i < num_vectors; i++){
    // //     Vector_pack(vectors[i], arrays[i]);
    // // }

    // char buffer[1000];
    // array_to_string(buffer, arrays[1], dim);
    // printf("array:     %s\n", buffer);

    // printf("Packing    ");
    // Vector_pack(vectors[1], arrays[1]);
    // printf("\n");
    // for(uint i = 0; i < dim; i++){
    //     printf("%d ", Vector_getEntry(vectors[1], i));
    // }
    // printf("\n");
    // vectorToString(buffer, vectors[1]);
    // printf("vector:    %s\n", buffer);
    // Vector_unpack(arrays[1], vectors[1]);
    // array_to_string(buffer, arrays[1], dim);
    // printf("unpacked:  %s\n", buffer);

    // return 0;
    // Vector_slice(vectors[1], vectors[0], 2, 4);


    // Vector *v = Vector_construct(&VectorGenericInterface, p, 4 - 2, vectors[1]->offset);
    // Vector_pack(v, arrays[3]);
    // vectorToString(buffer, v);
    // printf("v:               %s\n", buffer);

    // vectorToString(buffer, vectors[0]);
    // printf("vector: %s\n", buffer);
    // vectorToString(buffer, vectors[1]);
    // printf("slice:           %s\n", buffer);
    
    // addVectorsGeneric(vectors[1], v, 1);
    // vectorToString(buffer, vectors[1]);
    // printf("sum:             %s\n", buffer);
    // vectorToString(buffer, vectors[0]);
    // printf("vector: %s\n", buffer);    

    // printf("argc : %d\n", argc);
    // if(argc > 1){
    //     printf("Testing addVectors\n");
    //     for(uint repeats = 0; repeats < num_repeats; repeats ++){
    //         for(uint i = 0; i < num_vectors; i++){
    //             for(uint j = 0; j < num_vectors; j++){
    //                 if(i==j){
    //                     continue;
    //                 }
    //                 addVectorsGeneric(p, vectors[i], vectors[j], scale);
    //             }
    //         }
    //     }
    // } else {
    //     printf("Comparison\n");
    //     for(uint repeats = 0; repeats < num_repeats; repeats ++){
    //         for(uint i = 0; i < num_vectors; i++){
    //             for(uint j = 0; j < num_vectors; j++){
    //                 if(i==j){
    //                     continue;
    //                 }
    //                 for(uint k = 0; k < dim; k++){
    //                     arrays[i][k] += arrays[j][k] * scale;
    //                     arrays[i][k] = arrays[i][k] % p;
    //                 }
    //             }
    //         }   
    //     }
    // }

    // array_to_string(buffer, arrays[0], 10);
    // printf("%s\n\n", buffer);
    // vectorToString(buffer, p, vectors[0]);
    // printf("%s\n\n", buffer);
    
    // #define ROWS 5
    // #define COLS 10
    // uint p = 5;

    #define ROWS 4
    #define COLS 9
    uint p = 2;
    initializePrime(p);
    Matrix *M = constructMatrix(&VectorGenericInterface, p, ROWS, COLS);
    uint array[ROWS][COLS] =     {
        {1, 1, 0, 0, 0,   1, 0, 0, 0},
        {0, 0, 0, 0, 0,   0, 1, 0, 0},
        {0, 0, 0, 1, 0,   0, 0, 1, 0},
        {1, 0, 0, 0, 0,   0, 0, 0, 1}
    };
    // uint array[ROWS][COLS] = {
    //     {1, 0, 0, 4, 0, 2, 1, 3, 0, 0},
    //     {4, 3, 2, 3, 4, 3, 4, 2, 4, 4},
    //     {3, 0, 1, 2, 3, 4, 4, 3, 1, 3},
    //     {0, 1, 3, 0, 2, 1, 3, 0, 2, 3},
    //     {4, 4, 3, 4, 2, 4, 0, 0, 0, 3}
    // };
    for(uint i = 0; i < ROWS; i++){
        Vector_pack(M->matrix[i], array[i]);
    }        
    printMatrix(M);  


    int pivots[COLS];
    rowReduce(M, pivots, COLS - ROWS, COLS - ROWS);
    printf("M: ");
    printMatrix(M);
    return 0;
}
*/