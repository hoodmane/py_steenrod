//
// Created by Hood on 5/22/2019.
//

/*
 * Vector implementation details:
 * For each prime $p$ we have the following fields:
 *      bitlength -- how long each component of the vector is in bits
 *      bitmask   -- mask out the bottom bitlength bits of a uint64 (so this is 2^bitlength - 1)
 *      entriesPer64Bits -- how many components we fit into a 64 bit integer.
 * 
 * A vector has a list of "limbs". Each limb is a 64 bit integer which stores entriesPer64Bits components. I took this limb
 * nomenclature from the GMP large integer library. I think it is reasonably evocative.
 * 
 * Public vector fields:
 *      dimension -- The dimension of the vector
 *      size      -- The size of the vector in bytes. equal to number of limbs * sizeof(uint64).
 *      offset    -- The number of empty fields at the beginning of the vector. Most often zero.
 *                   Vectors with nonzero offset come from slicing a longer vector, or if you are preparing to write to a specific slice.
 * 
 * Private fields:
 *      vectorImplementation -- A bundle of function pointers for the arithmetic. The arithmetic operations are different at 2 vs not 2
 *                              because at 2 we can just bit xor to add, whereas with other primes we have to do a bit more work.
 *      number_of_limbs      -- A convenience field. Equal to ceil((dimension + offset/bit_length)/entriesPer64Bits), also equal to size/sizeof(uint64)
 *      limbs                -- The actual backing for the vector.
 */ 
          

#include "FpVector.h"
#include "combinatorics.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int array_toString(char *buffer, uint *A, uint length){
    if(length == 0){
        return sprintf(buffer, "[]"); 
    }
    buffer[0] = '[';
    buffer[1] = '\0';
    int len = 1;
    for (uint i = 0; i < length; i++) {
        len += sprintf(buffer + len, "%d, ", A[i]);
    }
    len -=2; // get rid of trailing ", "
    len += sprintf(buffer + len, "]");
    return len;
}

void array_print(char *format_string, uint *A, uint length){
    char buffer[1000];
    array_toString(buffer, A, length);
    printf(format_string, buffer);
}


typedef struct VectorImplementation_s VectorImplementation;
VectorImplementation *getVectorImplementation(uint p);

typedef struct {
    uint dimension; // These have to match the definition of Vector in FpVector.h
    uint size;
    uint offset;
// Private fields:
    struct VectorImplementation_s *implementation; // Function pointers to arithmetic. 
    uint number_of_limbs;
    uint64 *limbs;
} VectorPrivate;

// Our arithmetic function pointers.
struct VectorImplementation_s {  
    uint p;  
    void (*addBasisElement)(Vector *target, uint idx, uint c);
    void (*addArray)(Vector *target, uint *source, uint c);
    void (*add)(Vector *target, Vector *source, uint c);
    void (*scale)(Vector *target, uint c);   
};

// Generated with Mathematica:
//     bitlengths = Prepend[#,1]&@ Ceiling[Log2[# (# - 1) + 1 &[Prime[Range[2, 54]]]]]
// But for 2 it should be 1.
uint bitlengths[MAX_PRIME_INDEX] = { 
     1, 3, 5, 6, 7, 8, 9, 9, 9, 10, 10, 11, 11, 11, 12, 12, 12, 12, 13,     
     13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15,    
     15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16 
};

uint getBitlength(uint p){
    return bitlengths[prime_to_index_map[p]];
}

// This is 2^bitlength - 1.
// Generated with Mathematica:
//     2^bitlengths-1
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

// This is floor(64/bitlength).
// Generated with Mathematica:
//      Floor[64/bitlengths]
uint entries_per_64_bits[MAX_PRIME_INDEX] = {
    64, 21, 12, 10, 9, 8, 7, 7, 7, 6, 6, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4,  
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

uint getEntriesPer64Bits(uint p){
    return entries_per_64_bits[prime_to_index_map[p]];   
}

char *modplookuptable[MAX_PRIME_INDEX] = {0};

// Called by initializePrime
void initializeModpLookupTable(uint p){
    uint p_times_p_minus_1 = p*(p-1) + 1;
    char *table = malloc((p_times_p_minus_1 + 1) * sizeof(uint));
    for(uint i = 0; i <= p_times_p_minus_1; i++){
        table[i] = i % p;
    }
    modplookuptable[prime_to_index_map[p]] = table;
}

// n must be in the range 0 <= n <= p * (p-1)
// This is probably stupid to do by lookup -- integer division is pretty fast.
// One idea is to do the lookup for 16 bits at a time. This would tighten up our
// inner loop so we only have to split each 64 bit limb into four parts, 
// and having a 65kb lookup table is an acceptable cost.
// If we don't do that, we should just use %p because this table is pointless.
// Maybe the lookup costs more than %p either way.
uint modPLookup(uint p, uint n){
    return modplookuptable[prime_to_index_map[p]][n];
}


typedef struct {
    uint limb;
    uint bit_index;
} LimbBitIndexPair;

LimbBitIndexPair *limbBitIndexLookupTable[MAX_PRIME_INDEX] = {0};

/**
 * Called by initializePrime
 * This table tells us which limb and which bitfield of that limb to look for a given index of
 * the vector in.
 */
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
    size += sizeof(VectorPrivate);
    return size; // Need extra space for vector fields.
}

/**
 * Pad the size of the vector to insure that later slicing it in a certain location
 * will have an offset of 0.
 */
uint Vector_getPaddedDimension(uint p, uint dimension, uint offset){
    uint entries_per_limb = getEntriesPer64Bits(p);
    uint bit_length = getBitlength(p);
    return ((dimension + offset/bit_length + entries_per_limb - 1)/entries_per_limb)*entries_per_limb;
}

/**
 * Since we use a lot of short lived vectors, we want to be able to stack allocate them.
 * This routine helps us do that. To stack allocate a vector:
 *      size_t container_size = Vector_getContainerSize(p);
 *      size_t total_size = container_size + Vector_getSize(p, dimension, offset);
 *      char memory[total_size];
 *      myVector = Vector_initialize(p, memory, dimension, offset);
 */
Vector *Vector_initialize(uint p, char **memory, uint dimension, uint offset){
    assert(offset < 64); // Offset >= 64 leads to undefined behavior!!
    VectorImplementation *vectImpl = getVectorImplementation(p);
    VectorPrivate *v = (VectorPrivate *) *memory;
    v->implementation = vectImpl;
    v->dimension = dimension;
    v->size = Vector_getSize(p, dimension, offset);
    size_t limb_size = v->size - sizeof(VectorPrivate);
    v->offset = offset;
    v->number_of_limbs = (limb_size + sizeof(uint64) - 1)/sizeof(uint64);
    v->limbs = dimension == 0 ? NULL : (uint64*)(*memory + sizeof(VectorPrivate));
    memset(v->limbs, 0, limb_size);
    *memory += v->size;
    return (Vector*)v;
}

/**
 * Construct a vector
 */
Vector *Vector_construct(uint p, uint dimension, uint offset){
    size_t container_size = sizeof(VectorPrivate);
    size_t total_size = container_size + Vector_getSize(p, dimension, offset);
    char *memory = malloc(total_size);
    Vector *result = Vector_initialize(p, &memory, dimension, offset);
    return result;
}

void Vector_free(Vector *vector){
    free(vector);
}

/**
 *  Zeroes a vector. This is lazy and just uses a memset, so if you've sliced it out of 
 *  a bigger vector, it's going to change some of the fields of the bigger vector.
 */
void Vector_setToZero(Vector *target){
    assert(target->offset == 0); // setToZero doesn't handle slices right now.
    VectorPrivate *t = (VectorPrivate*) target;
    memset(t->limbs, 0, t->size - sizeof(VectorPrivate));
}

/**
 * Assign one vector to another.
 * This is lazy and just uses a memcpy. Again, it will change fields outside the end
 * of a vector if it's sliced from a larger one. The asserts ensure that no problem
 * happens at the start, but we don't have a way to check at the end.
 */
void Vector_assign(Vector *target, Vector *source){
    assert(source->dimension == target->dimension);
    assert(source->offset == target->offset);
    assert(source->offset == 0); // Vector_assign doesn't handle slices right now.
    VectorPrivate *t = (VectorPrivate*) target;
    VectorPrivate *s = (VectorPrivate*) source;
    memcpy(t->limbs, s->limbs, s->size - sizeof(VectorPrivate));
}

/**
 * Unpacks the limb at index limb_idx into limb_array. limb_array should be of length 
 * entriesPer64Bits.
 * 
 * These unpackLimb and packLimb are building blocks for all of our other methods.
 * The messy boundary logic is to make sure they respect slices correctly.
 */
uint unpackLimb(uint *limb_array, VectorPrivate *vector, uint limb_idx){
    uint bit_length = getBitlength(vector->implementation->p);
    uint entries_per_64_bits = getEntriesPer64Bits(vector->implementation->p);
    uint bit_mask = getBitMask(vector->implementation->p);    
    uint bit_min = 0;
    uint bit_max = bit_length * entries_per_64_bits;    
    if(limb_idx == 0){
        bit_min = vector->offset;
    }
    if(limb_idx == vector->number_of_limbs - 1){
        // Calculates how many bits of the last field we need to use. But if it divides
        // perfectly, we want bit max equal to bit_length * entries_per_64_bits, so subtract and add 1.
        // to make the output in the range 1 -- bit_length * entries_per_64_bits.
        uint bits_needed_for_entire_vector = vector->offset + vector->dimension * bit_length;
        uint usable_bits_per_limb = bit_length * entries_per_64_bits;
        bit_max = 1 + ((-1 + bits_needed_for_entire_vector)%(usable_bits_per_limb));
    }

    uint64 limb_value = vector->limbs[limb_idx];
    uint idx = 0;
    for(uint j = bit_min; j < bit_max; j += bit_length){
        limb_array[idx] = (limb_value >> j) & bit_mask;
        idx++;
    }
    return idx;
}

uint packLimb(VectorPrivate *vector, uint *limb_array, uint limb_idx){
    uint bit_length = getBitlength(vector->implementation->p);
    uint entries_per_64_bits = getEntriesPer64Bits(vector->implementation->p);
    uint bit_min = 0;
    uint bit_max = bit_length * entries_per_64_bits; 
    // Set up and handle first and last limb   
    if(limb_idx == 0){
        bit_min = vector->offset;
    }
    if(limb_idx == vector->number_of_limbs - 1){
        uint bits_needed_for_entire_vector = vector->offset + vector->dimension * bit_length;
        uint usable_bits_per_limb = bit_length * entries_per_64_bits;
        bit_max = 1 + ((-1 + bits_needed_for_entire_vector)%(usable_bits_per_limb));
    }
    uint64 bit_mask = 0;
    if(bit_max - bit_min < 64){
        // Critical that we're bitshifting 1LL here. Otherwise, shifting more than 32 bits
        // is undefined behavior!!
        bit_mask = (1LL << (bit_max - bit_min)) - 1;
        bit_mask <<= bit_min;
        bit_mask = ~bit_mask;
    }
    // copy data in
    uint idx = 0;
    uint64 limb_value = vector->limbs[limb_idx] & bit_mask;
    for(uint j = bit_min; j < bit_max; j += bit_length){
        limb_value |= ((uint64) limb_array[idx]) << j;
        idx ++;
    }
    vector->limbs[limb_idx] = limb_value;
    return idx;
}

// Just pack each limb individually.
void Vector_pack(Vector *target, uint *source){
    assert(target != NULL);
    VectorPrivate * t = (VectorPrivate*) target;
    for(uint limb = 0; limb < t->number_of_limbs; limb++){
        source += packLimb(t, source, limb);
    }
}

void Vector_unpack(uint * target, Vector * source){ 
    assert(target != NULL);
    VectorPrivate * s = (VectorPrivate*) source;    
    for(uint i = 0; i < s->number_of_limbs; i++){
        target += unpackLimb(target, s, i);
    }
}

// Find the limb index is located in with getLimbBitIndexPair.
// Bit shift down so the entry is the least significant part of the limb
// and mask it out.
uint Vector_getEntry(Vector *vector, uint index){
    assert(vector != NULL);
    assert(index < vector->dimension);
    VectorPrivate *v = (VectorPrivate*) vector;    
    uint bit_length = getBitlength(v->implementation->p);    
    uint64 bit_mask = getBitMask(v->implementation->p);
    LimbBitIndexPair limb_index = getLimbBitIndexPair(v->implementation->p, index + v->offset/bit_length);
    uint64 result = v->limbs[limb_index.limb];
    result >>= limb_index.bit_index;
    result &= bit_mask;
    return result;
}

void Vector_setEntry(Vector *vector, uint index, uint value){
    assert(vector != NULL);
    assert(index < vector->dimension);
    VectorPrivate *v = (VectorPrivate*) vector;    
    uint bit_length = getBitlength(v->implementation->p);    
    uint64 bit_mask = getBitMask(v->implementation->p);
    LimbBitIndexPair limb_index = getLimbBitIndexPair(v->implementation->p, index + v->offset/bit_length);
    uint64 *result = &(v->limbs[limb_index.limb]);
    *result &= ~(bit_mask << limb_index.bit_index);
    *result |= (((uint64)value) << limb_index.bit_index);
}

void Vector_slice(Vector *result, Vector *source, uint start, uint end){
    assert(start <= end);
    assert(end <= source->dimension);
    VectorPrivate *r = (VectorPrivate*) result;
    VectorPrivate *s = (VectorPrivate*) source;
    // Use the implementation field to detect erroneously declared slice targets.
    assert(r->implementation == s->implementation);    
    r->dimension = end - start;
    if(start == end){
        r->size = sizeof(VectorPrivate);
        r->number_of_limbs = 0;
        r->limbs = NULL;
        return;    
    }
    uint bit_length = getBitlength(s->implementation->p);
    LimbBitIndexPair limb_index = getLimbBitIndexPair(r->implementation->p, start + source->offset/bit_length);
    r->offset = limb_index.bit_index;
    r->size = Vector_getSize(r->implementation->p, r->dimension, r->offset);
    size_t limb_size = r->size - sizeof(VectorPrivate);
    r->number_of_limbs = limb_size/sizeof(uint64);
    r->limbs = s->limbs + limb_index.limb;
}

VectorIterator Vector_getIterator(Vector *vector){
    VectorIterator result;
    VectorPrivate *v = (VectorPrivate*)vector;
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
    result.value = ((v->limbs[result.limb_index]) >> result.bit_index) & bit_mask;
    return result;
}

VectorIterator Vector_stepIterator(VectorIterator it){
    it.index ++;
    VectorPrivate *v = (VectorPrivate*)it.vector;
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
    it.value = ((v->limbs[it.limb_index]) >> it.bit_index) & bit_mask;
    return it;
}

// For the arithmetic, at 2 we xor whereas at odd primes we have to do a bit more work.
// Odd prime vector arithmetic
// We've chosen our packing so that we have enough space to fit p*(p-1) in each limb.
// So we do the arithmetic in place. Then we have to reduce each slot by p.
// We do this by pulling the slots out and reducing each one separately mod p via table lookup.

void VectorGeneric_addBasisElement(Vector *target, uint index, uint coeff){
    assert(index < target->dimension);
    VectorPrivate *t = (VectorPrivate*) target;    
    uint64 bit_mask = getBitMask(t->implementation->p);
    uint bit_length = getBitlength(t->implementation->p);    
    LimbBitIndexPair limb_index = getLimbBitIndexPair(t->implementation->p, index + target->offset/bit_length);
    uint64 *result = &(t->limbs[limb_index.limb]);
    uint64 new_entry = *result >> limb_index.bit_index;
    new_entry &= bit_mask;
    new_entry += coeff;
    new_entry = modPLookup(t->implementation->p, new_entry);
    *result &= ~(bit_mask << limb_index.bit_index);
    *result |= (new_entry << limb_index.bit_index);
}

void VectorGeneric_addArray(Vector *target, uint *source, uint c){
    VectorPrivate *t = (VectorPrivate*) target;
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
    VectorPrivate *t = (VectorPrivate*) target;
    VectorPrivate *s = (VectorPrivate*) source;    
    uint entries[getEntriesPer64Bits(t->implementation->p)];    
    for(uint i = 0; i < s->number_of_limbs; i++){
        t->limbs[i] = t->limbs[i] + coeff * s->limbs[i];
        uint limb_length = unpackLimb(entries, t, i);
        for(uint j = 0; j < limb_length; j++){
            entries[j] = modPLookup(t->implementation->p, entries[j]);
        }
        packLimb(t, entries, i);
    }
}

void VectorGeneric_scale(Vector *target, uint coeff){
    // assert(coeff != 0);
    VectorPrivate *t = (VectorPrivate*) target;   
    uint entries_per_64_bits = getEntriesPer64Bits(t->implementation->p);    
    uint entries[entries_per_64_bits];
    for(uint i = 0; i < t->number_of_limbs; i++){
        t->limbs[i] = coeff * t->limbs[i];
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
    VectorPrivate *t = (VectorPrivate*) target;    
    uint bit_length = getBitlength(t->implementation->p);
    LimbBitIndexPair limb_index = getLimbBitIndexPair(t->implementation->p, index + target->offset/bit_length);
    uint64 *result = &(t->limbs[limb_index.limb]);
    *result ^= ((uint64)coeff << limb_index.bit_index);
}

void Vector2_addArray(Vector *target, uint *source, uint c __attribute__((unused))){
    VectorPrivate *t = (VectorPrivate*) target;
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
        t->limbs[limb] ^= source_limb;
    }
}

void Vector2_add(Vector *target, Vector *source, uint coeff){
    assert(target->dimension == source->dimension);
    assert(target->offset == source->offset);
    VectorPrivate *t = (VectorPrivate*) target;
    VectorPrivate *s = (VectorPrivate*) source;
    assert(t->implementation == s->implementation);
    uint bit_length = getBitlength(t->implementation->p); 
    uint64 *target_ptr = t->limbs;
    uint64 *source_ptr = s->limbs;
    if(t->number_of_limbs == 0){
        return;
    }
    uint64 source_limb;
    if(t->number_of_limbs == 1){
        source_limb = *source_ptr;
        uint64 bit_mask = -1;
        uint bit_min = source->offset;
        uint bit_max = (source->offset + source->dimension * bit_length) % 64;
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
    uint bit_max = (source->offset + source->dimension*bit_length) % 64;
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
    VectorPrivate *t = (VectorPrivate*) target;   
    uint64 *target_ptr = t->limbs;
    for(uint i = 0; i < t->number_of_limbs; i++){
        *target_ptr *= coeff;
        target_ptr++;
    }
}

uint64 Vector3_reduceLimb(uint64 limb);
// Vector3_add -- same as VectorGeneric_add except for the contents of Vector3_reduceLimb
void Vector3_add(Vector *target, Vector *source, uint coeff){
    assert(source->dimension == target->dimension);
    assert(target->offset == source->offset);
    VectorPrivate *t = (VectorPrivate*) target;
    VectorPrivate *s = (VectorPrivate*) source;       
    for(uint i = 0; i < s->number_of_limbs; i++){
        t->limbs[i] += coeff * s->limbs[i];
        t->limbs[i] = Vector3_reduceLimb(t->limbs[i]);
    }
}

// Vector3_scale -- same as VectorGeneric_scale except for the contents of Vector3_reduceLimb
void Vector3_scale(Vector *target, uint coeff){
    // assert(coeff != 0);
    VectorPrivate *t = (VectorPrivate*) target;   
    for(uint i = 0; i < t->number_of_limbs; i++){
        t->limbs[i] *= coeff;
        t->limbs[i] = Vector3_reduceLimb(t->limbs[i]);
    }
}

uint64 Vector3_reduceLimb(uint64 limb){
    uint64 top_bit_set_in_each_field = 0x4924924924924924LL;
    limb = ((limb & top_bit_set_in_each_field) >> 2) + (limb & (~top_bit_set_in_each_field));
    uint64 limb_3s = limb & (limb >> 1);
    limb_3s |= limb_3s << 1;
    limb ^= limb_3s;
    return limb;
}

uint64 Vector5_reduceLimb(uint64 limb);
// Vector3_add -- same as VectorGeneric_add except for the contents of Vector3_reduceLimb
void Vector5_add(Vector *target, Vector *source, uint coeff){
    assert(source->dimension == target->dimension);
    assert(target->offset == source->offset);
    VectorPrivate *t = (VectorPrivate*) target;
    VectorPrivate *s = (VectorPrivate*) source;       
    for(uint i = 0; i < s->number_of_limbs; i++){
        t->limbs[i] += coeff * s->limbs[i];
        t->limbs[i] = Vector5_reduceLimb(t->limbs[i]);
    }
}

// Vector3_scale -- same as VectorGeneric_scale except for the contents of Vector3_reduceLimb
void Vector5_scale(Vector *target, uint coeff){
    // assert(coeff != 0);
    VectorPrivate *t = (VectorPrivate*) target;   
    for(uint i = 0; i < t->number_of_limbs; i++){
        t->limbs[i] *= coeff;
        t->limbs[i] = Vector5_reduceLimb(t->limbs[i]);
    }
}

uint64 Vector5_reduceLimb(uint64 limb){
    uint64 bottom_bit = 0x84210842108421LL;
    uint64 bottom_two_bits = bottom_bit | (bottom_bit << 1);
    uint64 bottom_three_bits = bottom_bit | (bottom_two_bits << 1);
    uint64 a = (limb >> 2) & bottom_three_bits;
    uint64 b = limb & bottom_two_bits;
    uint64 m = (bottom_bit << 3) - a + b;
    uint64 c = (m >> 3) & bottom_bit;
    c |= c << 1;
    uint64 d = m & bottom_three_bits;
    return d + c - bottom_two_bits;



    // uint64 input_limb = limb;
    // // unset top two bits.
    // limb &= ~((bottom_bit_set_in_each_field << 4) | (bottom_bit_set_in_each_field << 3));
    // // 16 --> 1
    // limb += (input_limb >> 4) & bottom_bit_set_in_each_field;
    // // 8 --> 3
    // uint64 w = (input_limb >> 3) & bottom_bit_set_in_each_field;
    // limb += w | (w << 1); // temp | (temp << 1) is 3*temp.
    // // Now we've unset the two leading bits. Need to reduce the rest mod 5.
    // // Maximum possible value at this point is 00111 = 7, so we have to subtract 5 zero or one time.
    // // To figure it out, we implement "<5" with bitops. We xor with 100 and find the leading bit of the result.
    // // If the leading bit is in the 1's or 2's place, the number was 5/6/7 and we need to subtract five.
    // // So with only leading bit set, or with self shifted by 1. Now the 1's bit is set if the leading bit was the 1's or 2's place.
    // // Now the 1's  bit is set if and only if the field contained a 5/6/7. We or with ourself shifted by 2 to get 5 or 0.
    // w = limb ^ (bottom_bit_set_in_each_field << 2); // xor with 100
    // w |= w>>1;
    // w &= ~(w>>1); // Only leading bit set
    // w |= w>>1; // If 2's bit is set, make sure 1's bit is set too
    // w &= bottom_bit_set_in_each_field; // is 1's bit set?
    // w |= w << 2; // Now w contains a 5 or a 0 depending on whether the field is >= 5.
    // limb -= w;
    // return limb;
}



VectorImplementation VectorGenericImplementation = {
    0,    
    VectorGeneric_addBasisElement, VectorGeneric_addArray, VectorGeneric_add, VectorGeneric_scale,
};

VectorImplementation Vector2Implementation = {
    0,
    Vector2_addBasisElement, Vector2_addArray, Vector2_add, Vector2_scale,
};

VectorImplementation Vector3Implementation = {
    0,
    VectorGeneric_addBasisElement, VectorGeneric_addArray, Vector3_add, Vector3_scale,
};

VectorImplementation Vector5Implementation = {
    0,
    VectorGeneric_addBasisElement, VectorGeneric_addArray, Vector5_add, Vector5_scale,
};

// The generic methods depend on the implementation, so we look it up on the target.
void Vector_addBasisElement(Vector *target, uint idx, uint c){
    assert(c < ((VectorPrivate*)target)->implementation->p);
    ((VectorPrivate*)target)->implementation->addBasisElement(target, idx, c);
}

void Vector_addArray(Vector *target, uint *source, uint c){
    assert(c < ((VectorPrivate*)target)->implementation->p);
    ((VectorPrivate*)target)->implementation->addArray(target, source, c);
}

void Vector_add(Vector *target, Vector *source, uint c){
    assert(c < ((VectorPrivate*)target)->implementation->p);
    ((VectorPrivate*)target)->implementation->add(target, source, c);
}

void Vector_scale(Vector *target, uint c){
    assert(c < ((VectorPrivate*)target)->implementation->p);
    ((VectorPrivate*)target)->implementation->scale(target, c);
}


VectorImplementation vectorImplementationTable[MAX_PRIME_INDEX];
// Called by initializePrime
void initializeVectorImplementation(uint p){
    if(p == 2){
        vectorImplementationTable[prime_to_index_map[p]] = Vector2Implementation;
    } else if(p==3) {
        vectorImplementationTable[prime_to_index_map[p]] = Vector3Implementation;
    } else if(p==5) {
        vectorImplementationTable[prime_to_index_map[p]] = Vector5Implementation;
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

void Vector_print(char *fmt_string, Vector *v){
    char buffer[10000];
    Vector_toString(buffer, v);
    printf("%s\n", buffer);
}

void Vector_serialize(char **buffer, Vector *v){
    VectorPrivate * vect = (VectorPrivate*) v;
    size_t limbs_size = v->size - sizeof(VectorPrivate);
    memcpy(*buffer, v, sizeof(VectorPrivate));
    memcpy(*buffer + sizeof(VectorPrivate), vect->limbs, limbs_size);
    ((VectorPrivate*)*buffer)->limbs = NULL;
    *buffer += v->size;
}

Vector *Vector_deserialize(char **source_buffer){
    VectorPrivate *v = (VectorPrivate*)*source_buffer;
    v->limbs = (uint64*)(*source_buffer + sizeof(VectorPrivate));
    *source_buffer += v->size;
    return (Vector*)v;
}

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
    matrix->matrix = vector_ptr;
    for(uint row = 0; row < rows; row++){
        *vector_ptr = Vector_initialize(p, &container_ptr, columns, 0);
        vector_ptr ++; 
    }
    return matrix;
}

void Matrix_serialize(char **buffer, Matrix *M){
    size_t size = Matrix_getSize(M->p, M->rows, M->columns);
    memcpy(*buffer, M, size);
    buffer += size;
}

Matrix *Matrix_deserialize(char **buffer){
    Matrix *matrix = (Matrix*)*buffer;
    *buffer += sizeof(Matrix);
    Vector **vector_ptr = (Vector**)buffer;
    matrix->matrix = vector_ptr;
    uint rows = matrix->rows;    
    *buffer += rows * sizeof(Vector*);
    for(uint row = 0; row < rows; row++){
        *vector_ptr = Vector_deserialize(buffer);
        vector_ptr ++;
    }
    return matrix;
}

Matrix *Matrix_construct(uint p,  uint rows, uint columns)  {
    char *M = malloc(Matrix_getSize(p, rows, columns));
    // printf("columns: %d, rows: %d\n", columns, rows);
    return Matrix_initialize(M, p, rows, columns);
}

void Matrix_free(Matrix *M){
    free(M);
}


uint Matrix_getSliceSize(uint p, uint rows){
    return sizeof(Matrix) + rows*(sizeof(Vector*) + sizeof(VectorPrivate));
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
    Vector **matrix_ptr = (Vector**)result->matrix; 
    char *vector_ptr = (char*)(matrix_ptr + num_rows);
    Vector *initialized_vector_ptr;
    for(uint i = 0; i < num_rows; i++){
        initialized_vector_ptr = Vector_initialize(M->p, &vector_ptr, 0, 0);
        *matrix_ptr = initialized_vector_ptr;
        Vector_slice(initialized_vector_ptr, M->matrix[i], column_min, column_max);
        matrix_ptr ++;
    }
    assert(matrix_ptr == (Vector **)(result->matrix + num_rows));
    assert(vector_ptr == (char*)matrix_ptr + num_rows * sizeof(VectorPrivate));
    return result;
}

uint Matrix_toString(char *buffer, Matrix *M){
    int len = 0;
    len += sprintf(buffer + len, "    [\n");
    for(uint i = 0; i < M->rows; i++){
        len += sprintf(buffer + len, "        ");
        len += Vector_toString(buffer + len, M->matrix[i]);
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
        char slice_memory[sizeof(VectorPrivate)];
        char *slice_ptr = slice_memory;
        Vector *slice = Vector_initialize(M->p, &slice_ptr, 0, 0);
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



Kernel *Kernel_construct(uint p, uint rows, uint columns){
    assert(columns < MAX_DIMENSION);
    Kernel *k = malloc(
        sizeof(Kernel) 
        + columns * sizeof(uint)
        + Matrix_getSize(p, rows, columns)
    );
    k->column_to_pivot_row = (int*)(k + 1);
    k->kernel = Matrix_initialize((char*)(k->column_to_pivot_row + columns), p, rows, columns);
    return k;
}

// size_t Kernel_serialize(char *buffer, Kernel *kernel){

// }

// Kernel *Kernel_deserialize(){

// }


void Kernel_free(Kernel *k){
    free(k);
}

/**
int main(int argc, char *argv[]){
    initializePrime(3);
    Vector *vect = Vector_construct(3, 10, 0);
    Vector *slice = Vector_construct(3, 0, 0);
    Vector_slice(slice, vect, 3, 7);
    uint target[4];
    printf("slice->dim: %d\n", slice->dimension);
    Vector_unpack(target, slice);
    array_print("target: %s\n", target, 4);

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


    // #define ROWS 4
    // #define COLS 9
    // uint p = 2;
    // initializePrime(p);
    // Matrix *M = constructMatrix(p, ROWS, COLS);
    // uint array[ROWS][COLS] =     {
    //     {1, 1, 0, 0, 0,   1, 0, 0, 0},
    //     {0, 0, 0, 0, 0,   0, 1, 0, 0},
    //     {0, 0, 0, 1, 0,   0, 0, 1, 0},
    //     {1, 0, 0, 0, 0,   0, 0, 0, 1}
    // };
    // // uint array[ROWS][COLS] = {
    // //     {1, 0, 0, 4, 0, 2, 1, 3, 0, 0},
    // //     {4, 3, 2, 3, 4, 3, 4, 2, 4, 4},
    // //     {3, 0, 1, 2, 3, 4, 4, 3, 1, 3},
    // //     {0, 1, 3, 0, 2, 1, 3, 0, 2, 3},
    // //     {4, 4, 3, 4, 2, 4, 0, 0, 0, 3}
    // // };
    // for(uint i = 0; i < ROWS; i++){
    //     Vector_pack(M->matrix[i], array[i]);
    // }        
    // printMatrix(M);  


    // int pivots[COLS];
    // rowReduce(M, pivots, COLS - ROWS, COLS - ROWS);
    // printf("M: ");
    // printMatrix(M);
    // return 0;
}
//*/