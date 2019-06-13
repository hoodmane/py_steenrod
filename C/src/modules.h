#ifndef CSTEENROD_MODULES_H
#define CSTEENROD_MODULES_H

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "algebra.h"

#define MODULE_TYPE_FREE 0
#define MODULE_TYPE_FINITE_DIMENSIONAL 1
#define MODULE_TYPE_FINITELY_PRESENTED 2

typedef struct Module {
    uint p;
    Algebra *algebra;    
    uint type;
    int min_degree;
    int max_degree;        
// Methods:
    bool (*computeBasis)(struct Module *this, int degree);
    uint (*getDimension)(struct Module *this, int degree);
    void (*actOnBasis)(struct Module *this, Vector *result, uint coeff, int op_degree, uint op_index, int mod_degree, uint mod_index);
} Module;

#define module_computeBasis(module, degree) ((module)->computeBasis)(module, degree)
#define module_getDimension(module, degree) ((module)->getDimension)(module, degree)
#define module_actOnBasis(module, result, coeff, op_deg, op, r_deg, r) ((module)->actOnBasis)(module, result, coeff, op_deg, op, r_deg, r)



typedef struct {
    Module module;
    int max_basis_degree;
    uint *graded_dimension;
    // This goes input_degree --> output_degree --> operation --> input_index --> Vector
    Vector ****actions;
} FiniteDimensionalModule;

FiniteDimensionalModule *FiniteDimensionalModule_construct(
    Algebra *algebra, int min_degree, int max_basis_degree, uint *graded_dimension);

void FiniteDimensionalModule_free(FiniteDimensionalModule *module);
void FiniteDimensionalModule_setAction(
    FiniteDimensionalModule *module,
    int operation_degree, uint operation_idx,
    int input_degree, uint input_idx,
    uint *output    
);

void FiniteDimensionalModule_setActionVector(
    FiniteDimensionalModule *module,
    int operation_degree, uint operation_idx,
    int input_degree, uint input_idx,
    Vector *output
);

Vector *FiniteDimensionalModule_getAction(
    FiniteDimensionalModule *module,
    int operation_degree, uint operation_idx,
    int input_degree, uint input_idx
);


bool FiniteDimensionalModule_computeBasis(Module *this, int degree);
uint FiniteDimensionalModule_getDimension(Module *this, int degree);
void FiniteDimensionalModule_actOnBasis(
    Module *this, Vector *result, uint coeff,  
    int op_degree, uint op_index, int mod_degree, uint mod_index);

typedef struct {
    Module module;
    uint *number_of_generators_in_degree;
} FreeModule;

typedef struct {
    int operation_degree;
    uint operation_index;
    int generator_degree;
    uint generator_index;
} FreeModuleOperationGeneratorPair;

FreeModule *FreeModule_construct(Algebra *algebra, int min_degree, int max_degree);
void FreeModule_free(FreeModule *module);

bool FreeModule_computeBasis(Module *this, int degree);
uint FreeModule_getNumberOfGensInDegree(FreeModule *this, int degree);
uint FreeModule_getDimension(Module *this, int degree);
void FreeModule_actOnBasis(Module *this, Vector *result, uint coeff, 
    int op_degree, uint op_index, int mod_degree, uint mod_idx);

void FreeModule_ConstructBlockOffsetTable(FreeModule *M, int degree);
uint FreeModule_operationGeneratorToIndex(FreeModule *this, int op_deg, uint op_idx, int gen_deg, uint gen_idx);
FreeModuleOperationGeneratorPair FreeModule_indexToOpGen(FreeModule *this, int degree, uint index);

typedef struct {
    int *column_to_pivot_row;
    Matrix *kernel;
} Kernel;

typedef struct {
    FreeModule *source;
    Module *target;
    Vector ***outputs; // degree --> input_idx --> output
    int max_degree;
    int max_computed_degree;
    Matrix **coimage_to_image_isomorphism;
    Kernel **kernel; // This is redundant with the next module's coimage_to_image_iso
} FreeModuleHomomorphism;
//void initializeFreeModuleHomomorphism(FreeModuleHomomorphism *f, )

FreeModuleHomomorphism *FreeModuleHomomorphism_construct(FreeModule *source, Module *target, int max_degree);
void FreeModuleHomomorphism_free(FreeModuleHomomorphism *f);
void FreeModuleHomomorphism_setOutput(FreeModuleHomomorphism *f, int gen_degree, uint gen_index, Vector *output);
void FreeModuleHomomorphism_applyToBasisElement(FreeModuleHomomorphism *f, Vector *result, uint coeff, int input_degree, uint input_index);

void FreeModuleHomomorphism_AllocateSpaceForNewGenerators(FreeModuleHomomorphism *f, int degree, uint num_gens);

void FreeModuleHomomorphism_getMatrix(FreeModuleHomomorphism *f, Matrix *result, int degree);

Kernel *Kernel_construct(uint p, uint rows, uint columns);
void Kernel_free(Kernel *k);

#endif //CSTEENROD_MODULES_H