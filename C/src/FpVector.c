//
// Created by Hood on 5/22/2019.
//

#include "FpVector.h"
#include "combinatorics.h"

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


typedef struct {
    uint64 dimension;
    uint64 size;
    uint64* vector;
} VectorStd;

typedef struct {
    uint block;
    uint index;
} VectorBlockIndexPair;

typedef struct {
    uint total_dimension; 
    uint total_size;
    uint number_of_blocks;
    uint * block_dimension;
    uint * block_sizes; 
    VectorBlockIndexPair * index_to_block;
} VectorBlockStructure;

// We're going to cast from Vector * to BlockVector * so first two fields really should be dimension and size.
// This is true in this case because the first two fields of the BlockStructure are dimension and size.
typedef struct {
    VectorBlockStructure blockStructure;
    Vector* blocks;
} BlockVector;

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

uint * modplookuptable[MAX_PRIME_INDEX] = {0};

void initializeModpLookupTable(uint p){
    uint p_times_p_minus_1 = p*(p-1);
    uint * table = malloc((p_times_p_minus_1 + 1) * sizeof(uint));
    for(uint i = 0; i <= p_times_p_minus_1; i++){
        table[i] = i % p;
    }
    modplookuptable[prime_to_index_map[p]] = table;
}

// n must be in the range 0 <= n <= p * (p-1)
uint modPLookup(uint p, uint n){
    return modplookuptable[prime_to_index_map[p]][n];
}

typedef struct {
    uint limb;
    uint bit_index;
} LimbBitIndexPair;

LimbBitIndexPair * limbBitIndexLookupTable[MAX_PRIME_INDEX] = {0};

void initializeLimbBitIndexLookupTable(uint p){
    uint p_idx = prime_to_index_map[p];
    uint entries_per_limb = getEntriesPer64Bits(p);
    uint bit_length = getBitlength(p);
    LimbBitIndexPair * table = malloc(MAX_DIMENSION * sizeof(LimbBitIndexPair));
    for(uint i = 0; i < MAX_DIMENSION; i++){
        table[i].limb = i/entries_per_limb;
        table[i].bit_index = (i % entries_per_limb) * bit_length;
    }
    limbBitIndexLookupTable[p_idx] = table;
}

LimbBitIndexPair getLimbBitIndexPair(uint p, uint idx){
    return limbBitIndexLookupTable[prime_to_index_map[p]][idx];
}

uint getVectorSize(uint p, uint dimension){
    uint size = (dimension == 0) ? 0 : (getLimbBitIndexPair(p, dimension - 1).limb + 1);
    return size + 3; // Need extra space for vector fields.
}

Vector* initializeVector(uint p, uint64 * memory, uint dimension){
    VectorStd * v = (VectorStd *) memory;
    v->dimension = dimension;
    v->size = getLimbBitIndexPair(p, dimension - 1).limb + 1;
    v->vector = dimension == 0 ? NULL : (uint64*)(v + 1);
    memset(v->vector, 0, v->size * sizeof(uint64));
    return (Vector*)v;
}

// There is no case distinction between Vector2 and VectorGeneric for the functions that just 
// get and set values.
Vector * constructVector(uint p, uint dimension){
    printf("Dimension: %d\n", dimension);
    uint size = dimension == 0 ? 0 : getLimbBitIndexPair(p, dimension - 1).limb + 1;
    VectorStd * result = malloc(
        sizeof(VectorStd) + size * sizeof(uint64)
    );
    result->dimension = dimension;
    result->size = size;
    result->vector = dimension == 0 ? NULL : (uint64*)(result + 1);
    memset(result->vector, 0, size * sizeof(uint64));
    return (Vector*)result;
}

void freeVector(Vector * vector){
    free(vector);
}

void assignVector(uint p, Vector * target, Vector * source){
    VectorStd * t = (VectorStd*) target;
    VectorStd * s = (VectorStd*) source;
    memcpy(t->vector, s->vector, t->size * sizeof(uint64));
}

void packVector(uint p, Vector * target, uint * source){
    VectorStd * t = (VectorStd*) target;
    uint bit_length = getBitlength(p);
    uint source_idx = 0;
    for(uint limb = 0; limb < target->size; limb++){
        uint64 limb_value = 0;
        for(uint bit_idx = 0; bit_idx < 64 - bit_length + 1 && source_idx < target->dimension; bit_idx += bit_length){
            limb_value |= ((uint64)source[source_idx]) << bit_idx;
            source_idx ++;
        }
        t->vector[limb] = limb_value;
    }
}

void unpackVector(uint p, uint * target, Vector * source){
    VectorStd * s = (VectorStd*) source;    
    uint bit_mask = getBitMask(p);
    uint bit_length = getBitlength(p);
    uint target_idx = 0;
    for(uint limb = 0; limb < s->size; limb++){
        uint64 result = s->vector[limb];
        for(uint bit_idx = 0; bit_idx < 64 - bit_length + 1 && target_idx < source->dimension; bit_idx += bit_length){
            target[target_idx] = result & bit_mask;
            target_idx++;
            result >>= bit_length;
        }
    }
}

uint getVectorEntry(uint p, Vector * vector, uint index){
    VectorStd * v = (VectorStd*) vector;    
    // assert index < vector->dimension
    uint64 bit_mask = getBitMask(p);
    LimbBitIndexPair limb_index = getLimbBitIndexPair(p, index);
    uint64 result = v->vector[limb_index.limb];
    result >>= limb_index.bit_index;
    result &= bit_mask;
    return result;
}

void setVectorEntry(uint p, Vector * vector, uint index, uint value){
    VectorStd * v = (VectorStd*) vector;    
    // assert index < vector->dimension
    uint64 bit_mask = getBitMask(p);
    LimbBitIndexPair limb_index = getLimbBitIndexPair(p, index);
    uint64 *result = &(v->vector[limb_index.limb]);
    *result &= ~(bit_mask << limb_index.bit_index);
    *result |= (((uint64)value) << limb_index.bit_index);
}

uint vectorToString(char * buffer, uint p, Vector * vector){
    VectorStd * v = (VectorStd*) vector;    
    uint bit_length = getBitlength(p);
    uint bit_mask = getBitMask(p);
    uint len = 0;
    len += sprintf(buffer + len, "[");
    uint target_idx = 0;
    for(uint limb = 0; limb < v->size; limb++){
        uint64 result = v->vector[limb];
        for(uint bit_idx = 0; bit_idx < 64 - bit_length + 1 && target_idx < v->dimension; bit_idx += bit_length){
            len += sprintf(buffer + len, "%lld, ", result & bit_mask);
            target_idx++;
            result >>= bit_length;
        }
    }
    len += sprintf(buffer + len, "]");
    return len;
}

// For the arithmetic, at 2 we xor whereas at odd primes we have to do a bit more work.

// Generic vector arithmetic
// We've chosen our packing so that we have enough space to fit p*(p-1) in each limb.
// So we do the arithmetic in place. Then we have to reduce each slot by p.
// We do this by pulling the slots out and reducing each one separately mod p via table lookup.

void addBasisElementToVectorGeneric(uint p, Vector * target, uint index, uint coeff){
    VectorStd * t = (VectorStd*) target;    
    uint bit_mask = getBitMask(p);
    LimbBitIndexPair limb_index = getLimbBitIndexPair(p, index);
    uint64 *result = &(t->vector[limb_index.limb]);
    uint64 new_entry = *result >> limb_index.bit_index;
    new_entry &= bit_mask;
    new_entry*= coeff;
    new_entry = modPLookup(p, new_entry);
    *result &= ~(bit_mask << limb_index.bit_index);
    *result |= (new_entry << limb_index.bit_index);
}

void addArrayGeneric(uint p, Vector * target, uint * source, uint c){
    VectorStd * t = (VectorStd*) target;
    uint bit_mask = getBitMask(p);
    uint bit_length = getBitlength(p);
    uint source_idx = 0;
    for(uint limb = 0; limb < target->size; limb++){
        uint64 target_limb = t->vector[limb];
        uint64 result = 0;
        for(uint j = 0; j < 64 - bit_length + 1 && source_idx < target->dimension; j += bit_length){
            uint64 entry = (target_limb & bit_mask);
            entry += c*source[source_idx];
            uint64 entry_mod_p = modPLookup(p, entry);
            result |= entry_mod_p << j;
            target_limb >>= bit_length;
            source_idx ++;
        }  
        t->vector[limb] = result;
    }
}

void addVectorsGeneric(uint p, Vector * target, Vector * source, uint coeff){
    VectorStd * t = (VectorStd*) target;
    VectorStd * s = (VectorStd*) source;    
    uint bit_length = getBitlength(p);
    uint bit_mask = getBitMask(p);
    uint64 * target_ptr = t->vector;
    uint64 * source_ptr = s->vector;
    for(uint i = 0; i < target->size; i++){
        uint64 total = *target_ptr + coeff * (*source_ptr);
        uint64 result = 0;
        for(uint j = 0; j < 64 - bit_length + 1; j += bit_length){
            uint64 entry_mod_p = modPLookup(p, (total & bit_mask));
            result |= entry_mod_p << j;
            total >>= bit_length;
        }
        *target_ptr = result;
        target_ptr++;
        source_ptr++;
    }
}

void scaleVectorGeneric(uint p, Vector * target, uint coeff){
    VectorStd * t = (VectorStd*) target;   
    uint bit_length = getBitlength(p);
    uint bit_mask = getBitMask(p);
    uint64 * target_ptr = t->vector;
    for(uint i = 0; i < target->size; i++){
        uint64 total = coeff * (*target_ptr);
        uint64 result = 0;
        for(uint j = 0; j < 64 - bit_length + 1; j += bit_length){
            uint64 entry_mod_p = modPLookup(p, (total & bit_mask));
            result |= entry_mod_p << j;
            total >>= bit_length;
        }
        *target_ptr = result;
        target_ptr++;
    }
}

// Vector2

void addBasisElementToVector2(uint p, Vector * target, uint index, uint coeff){
    VectorStd * t = (VectorStd*) target;    
    uint bit_mask = getBitMask(p);
    LimbBitIndexPair limb_index = getLimbBitIndexPair(p, index);
    uint64 *result = &(t->vector[limb_index.limb]);
    uint64 new_entry = *result >> limb_index.bit_index;
    new_entry &= bit_mask;
    new_entry *= coeff;
    *result &= ~(bit_mask << limb_index.bit_index);
    *result |= (new_entry << limb_index.bit_index);
}

void addArray2(uint p, Vector * target, uint * source, uint c){
    VectorStd * t = (VectorStd*) target;
    uint bit_length = getBitlength(p);
    uint source_idx = 0;
    for(uint limb = 0; limb < target->size; limb++){
        uint64 source_limb = 0;
        for(uint j = 0; j < 64 - bit_length + 1 && source_idx < target->dimension; j += bit_length){
            source_limb |= source[source_idx] << j;
            source_idx ++;
        }  
        t->vector[limb] ^= source_limb;
    }
}

void addVectors2(uint p, Vector * target, Vector * source, uint coeff){
    VectorStd * t = (VectorStd*) target;
    VectorStd * s = (VectorStd*) source;    
    uint64 * target_ptr = t->vector;
    uint64 * source_ptr = s->vector;
    for(uint i = 0; i < target->size; i++){
        *target_ptr ^= (coeff * (*source_ptr));
        target_ptr++;
        source_ptr++;
    }
}

void scaleVector2(uint p, Vector * target, uint coeff){
    VectorStd * t = (VectorStd*) target;   
    uint64 * target_ptr = t->vector;
    for(uint i = 0; i < target->size; i++){
        *target_ptr *= coeff;
        target_ptr++;
    }
}

VectorInterface VectorGenericInterface = {
    assignVector, packVector, unpackVector, vectorToString,
    getVectorEntry, setVectorEntry,
    addBasisElementToVectorGeneric, addArrayGeneric, addVectorsGeneric, scaleVectorGeneric,
    getVectorSize, initializeVector,
    constructVector, freeVector
};

VectorInterface Vector2Interface = {
    assignVector, packVector, unpackVector, vectorToString,
    getVectorEntry, setVectorEntry,
    addBasisElementToVector2, addArray2, addVectors2, scaleVector2,
    getVectorSize, initializeVector,
    constructVector, freeVector
};

void printVector(uint p, Vector * v){
    char buffer[200];
    vectorToString(buffer, p, v);
    printf("%s\n", buffer);
}

/*
int main(){
    initializePrime(5);
    int p = 5;
    uint v_arr[10] = {1, 2, 3, 4, 0, 4, 3, 2, 1, 0};
    uint w_arr[10] = {4, 3, 2, 1, 1, 2, 3, 4, 1, 1};
    Vector * v = constructVector(p, 10);
    packVector(p, v, v_arr);
    Vector * w = constructVector(p, 10);
    packVector(p, w, w_arr);    
    char buffer[1000];
    vectorToString(buffer, p, v);
    printf("%s\n", buffer);
    addVectorsGeneric(5, v, w, 2);
    scaleVectorGeneric(5, w, 2);
    vectorToString(buffer, p, v);
    printf("%s\n", buffer);    
    vectorToString(buffer, p, w);
    printf("%s\n", buffer); 

    initializePrime(2);
    uint v2_arr[10] = {1, 1, 0, 1, 1, 0, 0, 1, 1, 0};
    uint w2_arr[10] = {1, 0, 1, 1, 1, 0, 0, 0, 1, 1};

    Vector * v2 = constructVector(2, 10);
    packVector(2, v2, v2_arr);
    Vector * w2 = constructVector(2, 10);
    packVector(2, w2, w2_arr);    
    vectorToString(buffer, 2, v2);
    printf("%s\n", buffer);
    vectorToString(buffer, 2, w2);
    printf("%s\n", buffer);    
    addVectors2(2, v2, w2, 1);    
    vectorToString(buffer, 2, v2);
    printf("%s\n", buffer);    
    return 0;
}
/**/

/*

// Invariants for a block structure:
//    total_dimension = total(block_dimensions)
//    total_size = total(block_sizes)
//    block_sizes[i] = ceil(block_dimensions[i] / entries_per_64_bits)
//    index_to_block 
VectorBlockStructure * constructVectorBlockStructure(uint p, uint number_of_blocks, uint * block_dimensions){
    uint total_dimension = 0;
    for(uint i = 0; i < number_of_blocks; i++){
        total_dimension += block_dimensions[i];
    }

    VectorBlockStructure * result = malloc(
            sizeof(VectorBlockStructure)
            + 2 * number_of_blocks * sizeof(uint) // block_dimensions and block_sizes
            + total_dimension * sizeof(VectorBlockIndexPair) // index_to_block
    );
    result->total_dimension = total_dimension;
    result->total_size = 0; // Don't know yet
    result->number_of_blocks = number_of_blocks;
    result->block_dimensions = (uint*)(result + 1);
    result->block_sizes = result->block_dimensions + number_of_blocks;
    result->index_to_block = result->block_sizes + number_of_blocks;

    memcpy(result->block_dimensions, block_dimensions, number_of_blocks * sizeof(uint));

    uint entries_per_64_bits = getEntriesPer64Bits(p);

    uint i = 0;
    for(uint block = 0; block < number_of_blocks; block++){
        result->block_sizes[block] = (result->block_dimensions[block]  + (entries_per_64_bits - 1)) / entries_per_64_bits;
        result->total_size += result->block_sizes[block];
        for(uint j = 0; j < block_dimensions[block]; j++){
            result->index_to_block[i] = (VectorBlockIndexPair){block, j};
            i++;
        }
    }
    return result;
}

void freeVectorBlockStructure(VectorBlockStructure * block_structure){
    free(block_structure);
};

// blockStructure should either be:
//   produced by constructBlockStructure
//   satisfy num_blocks = 1 and have the total_dimension field set.
Vector * constructVector(uint p, VectorBlockStructure * blockStructure){
    uint total_size = blockStructure->total_size;
    uint num_blocks = blockStructure->number_of_blocks;
    uint* block_sizes = blockStructure->block_sizes;
    Vector * result = malloc(sizeof(Vector) + num_blocks * sizeof(uint64*) + total_size * sizeof(uint64));
    result->blockStructure = *blockStructure;
    result->vector_blocks = (uint64**)(result + 1);
    uint64* block_ptr = (uint64 *) (result->vector_blocks + num_blocks);
    for(uint i = 0; i < num_blocks; i++){
        result->vector_blocks[i] = block_ptr;
        block_ptr += blockStructure->block_sizes[i];
    }
    // assert(block_ptr == (uint64 *) (result->vector_blocks + num_blocks) + total_size);
    return result;
}

void freeVector(Vector * vector){
    free(vector);
}



void assignVector(uint p, Vector * target, Vector * source){
    memcpy(&target->vector_blocks[0][0], &source->vector_blocks[0][0], target->blockStructure.total_size * sizeof(uint64));
}

// source should be a single block of size block_idx
void assignVectorBlock(uint p, Vector * target, Vector * source, uint block_idx){
    memcpy(target->vector_blocks[block_idx], source->vector_blocks[0], target->blockStructure.block_sizes[block_idx]);
}

// Produces a single block Vector which points to block_idx from source.
// No copy is performed and writing to this vector writes into the target block.
// The memory that this points to is owned by the larger vector.
// This vector still has to be freed but it doesn't own the memory where its contents are stored.
Vector * constructExtractedVectorBlock(uint p, Vector * source, uint block_idx){
    Vector * v = malloc(sizeof(Vector));
    v->blockStructure.total_dimension = source->blockStructure.block_dimensions[block_idx];
    v->blockStructure.total_size = source->blockStructure.block_sizes[block_idx];
    v->blockStructure.number_of_blocks = 1;
    v->blockStructure.block_dimensions = &v->blockStructure.total_dimension;
    v->blockStructure.block_sizes = &v->blockStructure.total_size;
    v->vector_blocks = &source[]
}

uint getVectorEntry(uint p, Vector * vector, uint index){
    // assert index < vector->blockStructure.total_dimension
    uint bit_length = getBitlength(p);
    uint bit_mask = getBitMask(p);
    VectorBlockIndexPair block_idx = vector->blockStructure.index_to_block[index];
    uint result = (vector->vector_blocks[block_idx.block][block_idx.index]);
    result >>= bit_length * block_idx.index;
    result &= bit_mask;
    return result;
}

void setVectorEntry(uint p, Vector * vector, uint index, uint value){
    // 
    uint bit_length = getBitlength(p);
    uint bit_mask = getBitMask(p);
    VectorBlockIndexPair block_idx = vector->blockStructure.index_to_block[index];
    uint *result = &(vector->vector_blocks[block_idx.block][block_idx.index]);
    *result &= ~(bit_mask << bit_length * block_idx.index);
    *result |= (value << bit_length * block_idx.index);
}


void addBasisElementToVector(uint p, Vector * target, uint index, uint coeff){
    uint bit_length = getBitlength(p);
    uint bit_mask = getBitMask(p);
    VectorBlockIndexPair block_idx = target->blockStructure.index_to_block[index];
    uint * result = &(target->vector_blocks[block_idx.block][block_idx.index]);
    uint new_entry = *result >> bit_length * block_idx.index;
    new_entry &= bit_mask;
    new_entry*= coeff;
    new_entry = modPLookup(p, new_entry);
    *result &= ~(bit_mask << bit_length * block_idx.index);
    *result |= (new_entry << bit_length * block_idx.index);
}

void addVectors(uint p, Vector * target, Vector * source, long coeff){
    uint bit_length = getBitlength(p);
    uint bit_mask = getBitMask(p);
    uint entries_per_64 = getEntriesPer64Bits(p);
    for(uint i = 0; i < target->blockStructure.number_of_blocks; i++){
        uint64 * target_ptr = target->vector_blocks[i];
        uint64 * source_ptr = source->vector_blocks[i];
        for(uint j = 0; j < target->blockStructure.block_sizes[i]; j++){
            uint64 total = *target_ptr + coeff * (*source_ptr);
            uint64 result = 0;
            for(uint k = 0; k < entries_per_64; k++){
                uint64 entry_mod_p = modPLookup(p, (total & bit_mask));
                result |= entry_mod_p << (bit_length * k);
                total >>= bit_length;
            }
            *target_ptr = result;
            source_ptr++;
            target_ptr++;
        }
    }
}

void scaleVector(uint p, Vector * target, uint coeff){
    uint bit_length = getBitlength(p);
    uint bit_mask = getBitMask(p);
    uint entries_per_64 = getEntriesPer64Bits(p);
    for(uint i = 0; i < target->blockStructure.number_of_blocks; i++){
        uint64 * target_ptr = target->vector_blocks[i];
        for(uint j = 0; j < target->blockStructure.block_sizes[i]; j++){
            uint64 total = coeff *  (*target_ptr);
            uint64 result = 0;
            for(uint k = 0; k < entries_per_64; k++){
                uint64 entry_mod_p = modPLookup(p, (total & bit_mask));
                result |=  entry_mod_p << (bit_length * k);
                total >>= bit_length;
            }
            *target_ptr = result;
            target_ptr++;
        }
    }
}

long vector_to_string(char * buffer, uint p, Vector * v){
    uint bit_length = getBitlength(p);
    uint bit_mask = getBitMask(p);
    uint entries_per_64 = getEntriesPer64Bits(p);    
    long len = 0;
    len += sprintf(buffer + len, "[");
    for(uint i = 0; i < v->blockStructure.number_of_blocks; i++){
        for(uint j = 0; j < v->blockStructure.block_sizes[i]; j++){
            uint64 block = v->vector_blocks[i][j];
            for(uint k = 0; k < entries_per_64; k++){
                len += sprintf(buffer + len, "%d, ", block & bit_mask);
                block >>= bit_length;
            }
        }
    }
    len += sprintf(buffer + len, "]");
    return len;
}


int** allocate_matrix(uint rows, uint cols)  {
    int** M = calloc(1, rows * sizeof(int*) + rows * cols * sizeof(int));
    for(int row = 0; row < rows; row++){
        M[row] = (int*)(M + rows) + row * cols;
    }
    return M;
}

int matrix_to_string(char * buffer, Vector *M, uint rows, uint columns){
    int len = 0;
    len += sprintf(buffer + len, "    [\n");
    for(int i = 0; i < rows; i++){
        len += vector_to_string(buffer + len, M[i]);
    }
    len += sprintf(buffer + len, "    ]\n");
    return len;
}

void print_matrix(uint p, Vector *matrix, uint rows, uint columns){
    char * buffer[10000];
    matrix_to_string(buffer, matrix, rows, columns);
    printf(buffer);
}

void row_reduce(uint p, Vector * matrix, int * column_to_pivot_row, uint rows, uint columns){
    memset(column_to_pivot_row, -1, columns * sizeof(uint));
    uint pivot = 0;
    for(uint pivot_column = 0; pivot_column < columns; pivot_column++){
        // Search down column for a nonzero entry.
        int pivot_row;
        for(pivot_row = pivot; pivot_row < rows; pivot_row ++){
            if(GetVectorEntry(matrix[pivot_row], pivot_column) != 0){
                break;
            }
        }
        // No pivot in pivot_column.
        if(pivot_row == rows){
            continue;
        }
        // Record position of pivot.
        column_to_pivot_row[pivot_column] = pivot_row;

        print_matrix(matrix, rows, columns);
        // Otherwise pivot_row contains a row with a pivot in current column.
        // swap pivot row up.
        int * temp = matrix[pivot];
        matrix[pivot] = matrix[pivot_row];
        matrix[pivot_row] = temp;

        printf("row(%ld) <==> row(%ld)\n", pivot, pivot_row);
        print_matrix(matrix, rows, columns);
        // Divide pivot row by pivot entry
        int c = GetVectorEntry(matrix[pivot], pivot_column);
        int c_inv = inverse(p, c);
        scaleVector(matrix[pivot], c_inv);

        printf("row(%ld) *= %ld\n", pivot, c_inv);
        print_matrix(matrix, rows, columns);
        for(int row = 0; row < rows; row++){
            // Between pivot and pivot_row, we already checked that the pivot column is 0, so skip ahead a bit.
            if(row == pivot){
                row = pivot_row;
                continue;
            }
            int pivot_column_entry = GetVectorEntry(matrix[row], pivot_column);
            int row_op_coeff = (p - pivot_column_entry) % p;
            if(row_op_coeff == 0){
                continue;
            }
            // Do row operation
            addVectors(matrix[row], matrix[pivot], row_op_coeff)
            print_matrix(matrix, rows, columns);
        }
        pivot ++;
    }
    return;
}
*/