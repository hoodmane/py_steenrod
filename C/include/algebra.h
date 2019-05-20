//
// Created by Hood on 5/20/2019.
//

#ifndef CSTEENROD_ALGEBRA_H
#define CSTEENROD_ALGEBRA_H
#include <stdbool.h>

// Careful with these macros: could cause multiple evaluation of algebra / module.
#define algebra_compute_basis(algebra, degree) (*(algebra)->compute_basis)(algebra, degree)
#define algebra_get_basis_dimension(algebra, degree) (*(algebra)->get_basis_dimension)(algebra, degree)
#define algebra_multiply_basis_elements(algebra, r_deg, r, s_deg, s) (*(algebra)->multiply_basis_elements)(algebra, r_deg, r, s_deg, s)

#define module_compute_basis(module, degree) (*(module)->compute_basis)(module, degree)
#define module_get_basis_dimension(module, degree) (*(module)->get_basis_dimension)(module, degree)
#define module_act_on_basis(module, degree, op_deg, op, r_deg, r) (*(module)->act_on_basis)(module, op_deg, op, r_deg, r)

typedef struct {
    unsigned long p;
    unsigned long degree;
    unsigned long dimension;
    long* vector;
} Vector;

typedef struct Algebra {
    unsigned long p;
// Methods:
    bool (*compute_basis)(struct Algebra* this, unsigned long degree);
    unsigned long (*get_basis_dimension)(struct Algebra* this, unsigned long degree);
    int (*multiply_basis_elements)(struct Algebra* this, unsigned long r_degree, unsigned long r, unsigned long s_degree, unsigned long s);
} Algebra;



typedef struct Module {
    unsigned long p;
    Algebra * algebra;
// Methods:
    bool (*compute_basis)(struct Module* this, unsigned long degree);
    unsigned long (*get_basis_dimension)(struct Module* this, unsigned long degree);
    int (*act_on_basis)(struct Module* this, unsigned long op_degree, unsigned long op_index, unsigned long mod_degree, unsigned long mod_index);
} Module;

Vector * allocateVector(unsigned long p, unsigned long degree, unsigned long dimension);
void freeVector(Vector * elt);

void addBasisElementToVector(Vector * elt, unsigned long idx, long coeff);
void addVector(Vector * target, Vector * source);
void scaleVector(Vector *, long);
void assignVector(Vector * target, Vector * source);

#endif //CSTEENROD_ALGEBRA_H
