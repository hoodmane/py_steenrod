//
// Created by Hood on 5/20/2019.
//

#ifndef CSTEENROD_ALGEBRA_H
#define CSTEENROD_ALGEBRA_H
#include <stdbool.h>

#include "FpVector.h"

/*
typedef struct {
    uint number_of_blocks;
    uint * block_sizes;
    uint * index_to_block;
    uint vector_size;
} VectorBlockStructure;

typedef struct {
    VectorBlockStructure * blockStructure;
    uint dimension;
    uint* vector;
} Vector;

void computeVectorBlockStructureSize(VectorBlockStructure * blockStructure);

Vector * allocateVector(uint p, VectorBlockStructure blockStructure);
void freeVector(Vector * vector);

uint getVectorEntry(uint p, Vector * v, uint index);
void setVectorEntry(uint p, Vector * v, uint index, uint value);

void addBasisElementToVector(uint p, Vector * elt, uint idx, uint c);
void addVectors(uint p, Vector * target, Vector * source, uint c);

void scaleVector(uint p, Vector * v, uint c);
void assignVector(uint p, Vector * target, Vector * source);


*/

typedef struct Algebra {
    uint p;
    VectorInterface vectorInterface;
// Methods:
    bool (*compute_basis)(struct Algebra* this, uint degree);
    uint (*get_dimension)(struct Algebra* this, uint degree);
    void (*multiply_basis_elements)(struct Algebra* this, Vector *result, uint coeff, uint r_degree, uint r, uint s_degree, uint s);
} Algebra;

// Careful with these macros: could cause multiple evaluation of algebra / module.
#define algebra_compute_basis(algebra, degree) (*(algebra)->compute_basis)(algebra, degree)
#define algebra_get_dimension(algebra, degree) (*(algebra)->get_dimension)(algebra, degree)
#define algebra_multiply_basis_elements(algebra, result, coeff, r_deg, r, s_deg, s) (*(algebra)->multiply_basis_elements)(algebra, result, coeff, r_deg, r, s_deg, s)



typedef struct Module {
    uint p;
    Algebra * algebra;
// Methods:
    bool (*compute_basis)(struct Module* this, uint degree);
    uint (*get_dimension)(struct Module* this, uint degree);
    void (*act_on_basis)(struct Module* this, Vector *result, uint coeff, uint op_degree,uint op_index, uint mod_degree, uint mod_index);
} Module;

#define module_compute_basis(module, degree) (*(module)->compute_basis)(module, degree)
#define module_get_basis_dimension(module, degree) (*(module)->getdimension)(module, degree)
#define module_act_on_basis(module, result, coeff, op_deg, op, r_deg, r) (*(module)->act_on_basis)(module, result, coeff, op_deg, op, r_deg, r)



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

/*
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
    FreeModule * source;
    Module * target;
    Vector ** outputs;
    uint max_computed_degree;
    Vector ** coimage_to_image_isomorphism;
} FreeModuleHomomorphism;
//void initializeFreeModuleHomomorphism(FreeModuleHomomorphism * f, )

void FreeModuleHomomorphism_apply_to_basis_element(FreeModuleHomomorphism * f, Vector * result, uint coeff, uint input_degree, uint input_index);

typedef struct {
    uint dimension;
    uint * column_to_pivot_row;
    Vector * kernel;
} Kernel;


// Resolution datatype
// We're storing the augmented resolution, the module we're resolving goes resolution_modules index 0
// (sort of -- it's not a FreeModule so we make a fake FreeModule with the appropriate dimension etc)
// Thus the index to resolution_modules is homological_degree + 1
// The index to resolution_differentials is (homogical_degree_of_source + 1).
typedef struct {
    Algebra * algebra; // These first two fields are probably not needed.
    Module * module;
    FreeModule * resolution_modules; // The index into resolution_modules is homological_degree + 1.
    FreeModuleHomomorphism * resolution_differentials;// Each differential has source the module with the same index in resolution_modules
    uint * internal_degree_to_resolution_stage;       // Records how far we've resolved in each degree (homological_degree + 1)
    Kernel * last_kernel;                   // Records the kernel of the highest differential computed in each degree (find out which by indexing internal_degree_to_resolution_stage)
} Resolution;


*/

#endif //CSTEENROD_ALGEBRA_H
