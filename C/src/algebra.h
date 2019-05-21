//
// Created by Hood on 5/20/2019.
//

#ifndef CSTEENROD_ALGEBRA_H
#define CSTEENROD_ALGEBRA_H
#include <stdbool.h>


typedef struct {
    unsigned long p;
    unsigned long degree;
    unsigned long dimension;
    long* vector;
} Vector;

Vector * allocateVector(unsigned long p, unsigned long degree, unsigned long dimension);
void freeVector(Vector * elt);

void addBasisElementToVector(Vector * elt, unsigned long idx, long coeff);
void addVectors(Vector * target, Vector * source);
void addArrayToVector(Vector * target, long * array);
void addVectorToArray(long * target, Vector * source);

void scaleVector(Vector *, long);
void assignVector(Vector * target, Vector * source);




typedef struct Algebra {
    unsigned long p;
// Methods:
    bool (*compute_basis)(struct Algebra* this, unsigned long degree);
    unsigned long (*get_basis_dimension)(struct Algebra* this, unsigned long degree);
    void (*multiply_basis_elements)(struct Algebra* this, Vector *result, unsigned long r_degree, unsigned long r, unsigned long s_degree, unsigned long s);
} Algebra;

// Careful with these macros: could cause multiple evaluation of algebra / module.
#define algebra_compute_basis(algebra, degree) (*(algebra)->compute_basis)(algebra, degree)
#define algebra_get_basis_dimension(algebra, degree) (*(algebra)->get_basis_dimension)(algebra, degree)
#define algebra_multiply_basis_elements(algebra, result, r_deg, r, s_deg, s) (*(algebra)->multiply_basis_elements)(algebra, result, r_deg, r, s_deg, s)



typedef struct Module {
    unsigned long p;
    Algebra * algebra;
// Methods:
    bool (*compute_basis)(struct Module* this, unsigned long degree);
    unsigned long (*get_basis_dimension)(struct Module* this, unsigned long degree);
    void (*act_on_basis)(struct Module* this, Vector *result, unsigned long op_degree, unsigned long op_index, unsigned long mod_degree, unsigned long mod_index);
} Module;

#define module_compute_basis(module, degree) (*(module)->compute_basis)(module, degree)
#define module_get_basis_dimension(module, degree) (*(module)->get_basis_dimension)(module, degree)
#define module_act_on_basis(module, result, op_deg, op, r_deg, r) (*(module)->act_on_basis)(module, result, op_deg, op, r_deg, r)



typedef struct {
    Module module;
    unsigned long dimension;
    unsigned long max_degree;
    unsigned long * number_of_basis_elements_in_degree;
    // This goes input_degree --> output_degree --> operation --> input_index --> Vector
    Vector **** actions;
} FiniteDimensionalModule;

FiniteDimensionalModule * constructFiniteDimensionalModule(Algebra * algebra, unsigned long dimension, unsigned long * generator_degrees);

void freeFiniteDimensionalModule(FiniteDimensionalModule * module);
void addActionToFiniteDimensionalModule(FiniteDimensionalModule * module,
                                        unsigned long operation_degree, unsigned long operation_idx,
                                        unsigned long input_degree, unsigned long input_idx,
                                        unsigned long * output
);

void FiniteDimensionalModule_act_on_basis(Module * this, Vector *result, unsigned long op_degree, unsigned long op_index, unsigned long mod_degree, unsigned long mod_index);


typedef struct {
    Module module;
    unsigned long number_of_generators;
    unsigned long * number_of_generators_in_degree;
} FreeModule;

typedef struct {
    FreeModule * source;
    Module * target;
    long ** outputs;
    unsigned long max_kernel_degree;
    long ** kernel;
} FreeModuleHomomorphism;

void FreeModuleHomomorphism_apply_to_basis_element(FreeModuleHomomorphism * f, Vector * result, unsigned long input_degree, unsigned long input_index);

typedef struct {
    Algebra * algebra;
    Module * module;
    FreeModule * resolution_modules;
    FreeModuleHomomorphism * resolution_differentials;
} Resolution;




#endif //CSTEENROD_ALGEBRA_H
