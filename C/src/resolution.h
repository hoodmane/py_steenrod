#ifndef CSTEENROD_RESOLUTION_H
#define CSTEENROD_RESOLUTION_H

#include "Algebra.h"
#include "Module.h"
#include "FiniteDimensionalModule.h"
#include "FreeModule.h"
#include "FreeModuleHomomorphism.h"

// Resolution datatype
// We're storing the augmented resolution, the module we're resolving goes resolution_modules index 0
// (sort of -- it's not a FreeModule so we make a fake FreeModule with the appropriate dimension etc)
// Thus the index to resolution_modules is homological_degree + 1
// The index to resolution_differentials is (homogical_degree_of_source + 1).
typedef struct Resolution_s {
    Algebra *algebra;
    Module *module;
    void (*addClass)(uint hom_deg, int int_deg, char *cocycle_name);
    void (*addStructline)(
        uint source_hom_deg, int source_int_deg, uint source_idx, 
        uint target_hom_deg, int target_int_deg, uint target_idx
    );
    uint max_homological_degree;
    int min_degree;
    int max_degree;
    FreeModule **modules; // The index into resolution_modules is homological_degree + 1.
    FreeModuleHomomorphism **differentials;// Each differential has source the module with the same index in resolution_modules
    int *internal_degree_to_resolution_stage;       // Records how far we've resolved in each degree (homological_degree + 1)
} Resolution;

// Getters for javascript
FreeModuleHomomorphism *Resolution_getDifferential(Resolution *resolution, uint homological_degree);


Resolution *Resolution_construct(
    FiniteDimensionalModule *module,
    int max_degree,
    void (*addClass)(uint hom_deg, int int_deg, char *cocycle_name),
    void (*addStructline)(
        uint source_hom_deg, int source_int_deg, uint source_idx, 
        uint target_hom_deg, int target_int_deg, uint target_idx
    )    
);
void Resolution_free(Resolution *resolution);

void Resolution_step(Resolution *resolution, uint homological_degree, int degree);
void Resolution_resolveThroughDegree(Resolution *res, int degree);

void Resolution_serialize(Resolution *res);

#endif // CSTEENROD_RESOLUTION_H