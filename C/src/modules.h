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
    uint max_degree; 
    uint degree_shift;       
// Methods:
    bool (*computeBasis)(struct Module *this, uint degree);
    uint (*getDimension)(struct Module *this, uint degree);
    void (*actOnBasis)(struct Module *this, Vector *result, uint coeff, uint op_degree,uint op_index, uint mod_degree, uint mod_index);
} Module;

#define module_computeBasis(module, degree) ((module)->computeBasis)(module, degree)
#define module_getDimension(module, degree) ((module)->getDimension)(module, degree)
#define module_actOnBasis(module, result, coeff, op_deg, op, r_deg, r) ((module)->actOnBasis)(module, result, coeff, op_deg, op, r_deg, r)



typedef struct {
    Module module;
    uint max_basis_degree;
    uint *graded_dimension;
    uint degree_shift;
    // This goes input_degree --> output_degree --> operation --> input_index --> Vector
    Vector ****actions;
} FiniteDimensionalModule;

FiniteDimensionalModule *FiniteDimensionalModule_construct(Algebra *algebra, uint max_generator_degree, uint *graded_dimension);

void FiniteDimensionalModule_free(FiniteDimensionalModule *module);
void FiniteDimensionalModule_setAction(
    FiniteDimensionalModule *module,
    uint operation_degree, uint operation_idx,
    uint input_degree, uint input_idx,
    uint *output    
);

void FiniteDimensionalModule_setActionVector(
    FiniteDimensionalModule *module,
    uint operation_degree, uint operation_idx,
    uint input_degree, uint input_idx,
    Vector *output
);

bool FiniteDimensionalModule_computeBasis(Module *this, uint dimension);
uint FiniteDimensionalModule_getDimension(Module *this, uint degree);
void FiniteDimensionalModule_actOnBasis(
    Module *this, Vector *result, uint coeff,  uint op_degree, uint op_index, uint mod_degree, uint mod_index);

typedef struct {
    Module module;
    uint *number_of_generators_in_degree;
} FreeModule;

typedef struct {
    uint operation_degree;
    uint operation_index;
    uint generator_degree;
    uint generator_index;
} FreeModuleOperationGeneratorPair;

FreeModule *FreeModule_construct(Algebra *algebra, uint max_degree);
void FreeModule_free(FreeModule *module);

bool FreeModule_computeBasis(Module *this, uint degree);
uint FreeModule_getDimension(Module *this, uint degree);
void FreeModule_actOnBasis(Module *this, Vector *result, uint coeff, uint op_degree, uint op_index, uint mod_degree, uint mod_idx);

void FreeModule_ConstructBlockOffsetTable(FreeModule *M, uint degree);
uint FreeModule_operationGeneratorToIndex(FreeModule *this, uint op_deg, uint op_idx, uint gen_deg, uint gen_idx);
FreeModuleOperationGeneratorPair FreeModule_indexToOpGen(FreeModule *this, uint degree, uint index);

typedef struct {
    int *column_to_pivot_row;
    Matrix *kernel;
} Kernel;

typedef struct {
    FreeModule *source;
    Module *target;
    Vector ***outputs; // degree --> input_idx --> output
    uint max_degree;
    uint max_computed_degree;
    Matrix **coimage_to_image_isomorphism;
    Kernel **kernel; // This is redundant with the next module's coimage_to_image_iso
} FreeModuleHomomorphism;
//void initializeFreeModuleHomomorphism(FreeModuleHomomorphism *f, )

FreeModuleHomomorphism *FreeModuleHomomorphism_construct(FreeModule *source, Module *target, uint max_degree);
void FreeModuleHomomorphism_free(FreeModuleHomomorphism *f);
void FreeModuleHomomorphism_setOutput(FreeModuleHomomorphism *f, uint gen_degree, uint gen_index, Vector *output);
void FreeModuleHomomorphism_applyToBasisElement(FreeModuleHomomorphism *f, Vector *result, uint coeff, uint input_degree, uint input_index);

void FreeModuleHomomorphism_AllocateSpaceForNewGenerators(FreeModuleHomomorphism *f, uint degree, uint num_gens);

void FreeModuleHomomorphism_getMatrix(FreeModuleHomomorphism *f, Matrix *result, uint degree);

Kernel *Kernel_construct(uint p, uint rows, uint columns);
void Kernel_free(Kernel *k);

#endif //CSTEENROD_MODULES_H