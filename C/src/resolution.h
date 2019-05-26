#ifndef CSTEENROD_RESOLUTION_H
#define CSTEENROD_RESOLUTION_H

#include "algebra.h"
#include "modules.h"

// Resolution datatype
// We're storing the augmented resolution, the module we're resolving goes resolution_modules index 0
// (sort of -- it's not a FreeModule so we make a fake FreeModule with the appropriate dimension etc)
// Thus the index to resolution_modules is homological_degree + 1
// The index to resolution_differentials is (homogical_degree_of_source + 1).
typedef struct {
    Algebra * algebra;
    Module * module;
    FreeModule * resolution_modules; // The index into resolution_modules is homological_degree + 1.
    FreeModuleHomomorphism * resolution_differentials;// Each differential has source the module with the same index in resolution_modules
    uint * internal_degree_to_resolution_stage;       // Records how far we've resolved in each degree (homological_degree + 1)
} Resolution;

#endif // CSTEENROD_RESOLUTION_H