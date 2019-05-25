//
// Created by Hood on 5/22/2019.
//

#ifndef C_FPVECTOR_H
#define C_FPVECTOR_H

#define MAX_DIMENSION 10000

typedef unsigned long long uint64;
typedef long long int64;
typedef unsigned int uint;

void initializeModpLookupTable(uint p);
void initializeLimbBitIndexLookupTable(uint p);

uint modPLookup(uint p, uint n);

    // Vector * constructVector(uint p, uint dimension);
    // void freeVector(Vector * v);

    // void assignVector(uint p, Vector * target, Vector * source);

    // uint getVectorEntry(uint p, Vector * v, uint index);
    // void setVectorEntry(uint p, Vector * v, uint index, uint value);

    // void addBasisElementToVector(uint p, Vector * elt, uint idx, uint c);
    // void addVectors(uint p, Vector * target, Vector * source, long c);
    // void scaleVector(uint p, Vector * v, uint c);

    // uint vectorToString(char * buffer, uint p, Vector * v);

typedef struct  {
    uint64 dimension;
    uint64 size;
} Vector;

typedef struct {    
    void (*assign)(uint p, Vector * target, Vector * source);
    void (*pack)(uint p, Vector * target, uint * source);
    void (*unpack)(uint p, uint * target, Vector * source);
    uint (*toString)(char * buffer, uint p, Vector * v);

    uint (*getEntry)(uint p, Vector * v, uint index);
    void (*setEntry)(uint p, Vector * v, uint index, uint value);
    
    void (*addBasisElement)(uint p, Vector * target, uint idx, uint c);
    void (*addArray)(uint p, Vector * target, uint * source, uint c);
    void (*add)(uint p, Vector * target, Vector * source, uint c);
    void (*scale)(uint p, Vector * target, uint c);

    uint (*getSize)(uint p, uint dimension);
    Vector * (*initialize)(uint p, uint64 * memory, uint dimension);
    Vector * (*construct)(uint p, uint dimension);
    void (*free)(Vector * v);    
} VectorInterface;

typedef struct {
    void (*assign)(uint p, Vector * target, Vector * source);
    void (*pack)(uint p, Vector * target, uint * source);
    void (*unpack)(uint p, uint * target, Vector * source);
    uint (*toString)(char * buffer, uint p, Vector * v);

    uint (*getEntry)(uint p, Vector * v, uint index);
    void (*setEntry)(uint p, Vector * v, uint index, uint value);
    
    void (*addBasisElement)(uint p, Vector * target, uint idx, uint c);
    void (*add)(uint p, Vector * target, Vector * source, uint c);
    void (*scale)(uint p, Vector * target, uint c);
} PartialVectorInterface;

extern VectorInterface VectorGenericInterface;
extern VectorInterface Vector2Interface;
//extern VectorInterface BlockVectorInterface;

Vector * constructVector(uint p, uint dimension);
void freeVector(Vector * v);

void assignVector(uint p, Vector * target, Vector * source);
void packVector(uint p, Vector * target, uint * source);
void unpackVector(uint p, uint * target, Vector * source);
uint vectorToString(char * buffer, uint p, Vector * v);

uint getVectorEntry(uint p, Vector * v, uint index);
void setVectorEntry(uint p, Vector * v, uint index, uint value);

void addBasisElementToVectorGeneric(uint p, Vector * elt, uint idx, uint c);
void addVectorsGeneric(uint p, Vector * target, Vector * source, uint c);
void scaleVectorGeneric(uint p, Vector * v, uint c);

void addBasisElementToVector2(uint p, Vector * elt, uint idx, uint c);
void addVectors2(uint p, Vector * target, Vector * source, uint c);
void scaleVector2(uint p, Vector * v, uint c);


/*
typedef struct BlockVector;
typedef struct BlockStructure;
VectorBlockStructure * constructVectorBlockStructure(uint p, uint number_of_blocks, uint * block_dimensions);
void freeVectorBlockStructure(VectorBlockStructure * vectorBlockStructure);

Vector * constructBlockVector(uint p, VectorBlockStructure * blockStructure);
void freeBlockVector(Vector * vector);


void printVector(uint p, Vector * v);


int** allocate_matrix(uint rows, uint cols);

void row_reduce(uint p, Vector * matrix, int * column_to_pivot_row, uint rows, uint columns);
*/

#endif //C_FPVECTOR_H
