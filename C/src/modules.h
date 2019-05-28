#ifndef CSTEENROD_MODULES_H
#define CSTEENROD_MODULES_H

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "algebra.h"

typedef struct Module {
    uint p;
    Algebra *algebra;
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
    uint dimension;
    uint max_degree;
    uint *number_of_basis_elements_in_degree;
    // This goes input_degree --> output_degree --> operation --> input_index --> Vector
    Vector ****actions;
} FiniteDimensionalModule;

FiniteDimensionalModule * FiniteDimensionalModule_construct(Algebra *algebra, uint dimension, uint *generator_degrees);

void FiniteDimensionalModule_free(FiniteDimensionalModule *module);
void FiniteDimensionalModule_setAction(
    FiniteDimensionalModule *module,
    uint operation_degree, uint operation_idx,
    uint input_degree, uint input_idx,
    Vector *output
);

bool FiniteDimensionalModule_computeBasis(Module *this, uint dimension);
uint FiniteDimensionalModule_getDimension(Module* this, uint degree);
void FiniteDimensionalModule_actOnBasis(
    Module *this, Vector *result, uint coeff,  uint op_degree, uint op_index, uint mod_degree, uint mod_index);

typedef struct {
    Module module;
    uint max_generator_degree;
    uint max_degree;
    uint number_of_generators;
    uint *number_of_generators_in_degree;
} FreeModule;

FreeModule *FreeModule_construct(Algebra * algebra, uint max_generator_degree, uint max_degree);
void FreeModule_free(FreeModule * module);

void FreeModule_ConstructBlockOffsetTable(FreeModule * M, uint degree);
uint FreeModule_operationGeneratorToIndex(FreeModule *this, uint op_deg, uint op_idx, uint gen_deg, uint gen_idx);

bool FreeModule_computeBasis(Module* this, uint degree);
uint FreeModule_getDimension(Module* this, uint degree);
void FreeModule_actOnBasis(Module * this, Vector * result, uint coeff, uint op_degree, uint op_index, uint mod_degree, uint mod_idx);

typedef struct {
    uint *column_to_pivot_row;
    Matrix *kernel;
} Kernel;

typedef struct {
    FreeModule *source;
    Module *target;
    Vector ***outputs;
    uint max_computed_degree;
    Matrix **coimage_to_image_isomorphism;
    Kernel **kernel; // This is redundant with the next module's coimage_to_image_iso
} FreeModuleHomomorphism;
//void initializeFreeModuleHomomorphism(FreeModuleHomomorphism * f, )

FreeModuleHomomorphism *constructFreeModuleHomomorphism(FreeModule *source, Module *target);
void FreeModuleHomomorphism_addGenerator(FreeModuleHomomorphism *f, uint gen_degree, Vector *output);
void FreeModuleHomomorphism_applyToBasisElement(FreeModuleHomomorphism *f, Vector *result, uint coeff, uint input_degree, uint input_index);

void FreeModuleHomomorphism_AllocateSpaceForNewGenerators(FreeModuleHomomorphism *f, uint num_gens);

void FreeModuleHomomorphism_getMatrix(Matrix *result, FreeModuleHomomorphism *f, uint degree);
Kernel *Kernel_construct(VectorInterface *vectImpl, uint p, uint rows, uint columns);

#endif //CSTEENROD_MODULES_H