#ifndef CSTEENROD_MODULES_H
#define CSTEENROD_MODULES_H

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "algebra.h"

typedef struct Module {
    uint p;
    Algebra * algebra;
// Methods:
    bool (*compute_basis)(struct Module* this, uint degree);
    uint (*get_dimension)(struct Module* this, uint degree);
    void (*act_on_basis)(struct Module* this, Vector *result, uint coeff, uint op_degree,uint op_index, uint mod_degree, uint mod_index);
} Module;

#define module_compute_basis(module, degree) ((module)->compute_basis)(module, degree)
#define module_get_dimension(module, degree) ((module)->get_dimension)(module, degree)
#define module_act_on_basis(module, result, coeff, op_deg, op, r_deg, r) ((module)->act_on_basis)(module, result, coeff, op_deg, op, r_deg, r)



typedef struct {
    Module module;
    uint dimension;
    uint max_degree;
    uint * number_of_basis_elements_in_degree;
    // This goes input_degree --> output_degree --> operation --> input_index --> Vector
    Vector **** actions;
} FiniteDimensionalModule;

FiniteDimensionalModule * constructFiniteDimensionalModule(Algebra * algebra, uint dimension, uint * generator_degrees);

void freeFiniteDimensionalModule(FiniteDimensionalModule * module);
void addActionToFiniteDimensionalModule(
    FiniteDimensionalModule * module,
    uint operation_degree, uint operation_idx,
    uint input_degree, uint input_idx,
    uint * output
);

void FiniteDimensionalModule_act_on_basis(Module * this, Vector *result, uint coeff,  uint op_degree, uint op_index, uint mod_degree, uint mod_index);


// This should be a private struct I think.
typedef struct {
    uint operation_degree;
    uint operation_index;
    uint generator_degree;
    uint generator_index;
} FreeModuleOperationGeneratorPair;

typedef struct {
    Module module;
    uint max_generator_degree;
    uint max_degree;
    uint number_of_generators;
    uint * number_of_generators_in_degree;
// private fields
    uint computed_degree;
    FreeModuleOperationGeneratorPair ** basis_element_to_opgen_table;
    uint ** generator_to_index_table;
} FreeModule;




typedef struct {
    uint dimension;
    uint * column_to_pivot_row;
    Vector ** kernel;
} Kernel;

typedef struct {
    FreeModule * source;
    Module * target;
    Vector *** outputs;
    uint max_computed_degree;
    Vector *** coimage_to_image_isomorphism;
    Kernel ** kernel; // This is redundant with the next module's coimage_to_image_iso
} FreeModuleHomomorphism;
//void initializeFreeModuleHomomorphism(FreeModuleHomomorphism * f, )

FreeModuleHomomorphism * constructFreeModuleHomomorphism(FreeModule * source, Module * target);
void addGeneratorToFreeModuleHomomorphism(FreeModuleHomomorphism * f, uint gen_degree, Vector * output);
void FreeModuleHomomorphism_apply_to_basis_element(FreeModuleHomomorphism * f, Vector * result, uint coeff, uint input_degree, uint input_index);

void FreeModuleAllocateSpaceForNewGenerators(FreeModuleHomomorphism * f, uint num_gens);

void getHomomorphismMatrix(Vector ** result, FreeModuleHomomorphism * f, uint degree);
Kernel * constructKernel(VectorInterface * vectImpl, uint p, uint rows, uint columns);

#endif //CSTEENROD_MODULES_H