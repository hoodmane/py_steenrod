#ifndef CSTEENROD_FINITE_DIMENSIONAL_MODULE_H
#define CSTEENROD_FINITE_DIMENSIONAL_MODULE_H
#include "Module.h"

typedef struct {
    Module module;
    int max_basis_degree;
    uint *graded_dimension;
    // This goes input_degree --> output_degree --> operation --> input_index --> Vector
    Vector *****actions;
} FiniteDimensionalModule;

FiniteDimensionalModule *FiniteDimensionalModule_construct(
    Algebra *algebra, int min_degree, int max_degree, uint *graded_dimension);

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

#endif //CSTEENROD_FINITE_DIMENSIONAL_MODULE_H