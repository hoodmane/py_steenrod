#include "FreeModule.h"
#include "FreeModuleHomomorphism.h"
#include "Resolution.h"

typedef struct {
    Resolution *source;
    Resolution *target;
    FreeModuleHomomorphism **maps;
    uint max_homological_degree;
    uint homological_degree_shift;
    int internal_degree_shift;
    int max_degree;
    int *computed_degree;
} ResolutionHomomorphism;

ResolutionHomomorphism *ResolutionHomomorphism_construct(
    Resolution *source, Resolution *target, 
    uint homological_degree_shift, int internal_degree_shift
);

void ResolutionHomomorphism_setBaseMap(ResolutionHomomorphism *f, int input_degree, int input_index, Vector *output);
void ResolutionHomomorphism_baseMapReady(ResolutionHomomorphism *f, int degree);
void ResolutionHomomorphism_extend(ResolutionHomomorphism *f, uint source_homological_degree, int source_degree);

FreeModuleHomomorphism *ResolutionHomomorphism_getMap(ResolutionHomomorphism *f, uint homological_degree);