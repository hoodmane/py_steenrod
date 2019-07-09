#ifndef CSTEENROD_FREE_MODULE_HOMOMORPHISM_H
#define CSTEENROD_FREE_MODULE_HOMOMORPHISM_H
#include "FreeModule.h"

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
void FreeModuleHomomorphism_applyToGenerator(FreeModuleHomomorphism *f, Vector *result, uint coeff, int generator_degree, uint generator_index);
void FreeModuleHomomorphism_applyToBasisElement(FreeModuleHomomorphism *f, Vector *result, uint coeff, int input_degree, uint input_index);

void FreeModuleHomomorphism_AllocateSpaceForNewGenerators(FreeModuleHomomorphism *f, int degree, uint num_gens);

void FreeModuleHomomorphism_getMatrix(FreeModuleHomomorphism *f, Matrix *result, int degree);


// For Javascript:
FreeModule *FreeModuleHomomorphism_getSource(FreeModuleHomomorphism *f);
Module *FreeModuleHomomorphism_getTarget(FreeModuleHomomorphism *f);

#endif //CSTEENROD_FREE_MODULE_HOMOMORPHISM_H