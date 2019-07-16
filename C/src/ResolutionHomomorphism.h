#ifndef CSTEENROD_RESOLUTION_HOMOMORPHISM_H
#define CSTEENROD_RESOLUTION_HOMOMORPHISM_H

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

void ResolutionHomomorphism_extendBaseMap(ResolutionHomomorphism *f, int degree);
void ResolutionHomomorphism_setBaseMap(ResolutionHomomorphism *f, int input_degree, int input_index, Vector *output);
void ResolutionHomomorphism_baseMapReady(ResolutionHomomorphism *f, int degree);
void ResolutionHomomorphism_extend(ResolutionHomomorphism *f, uint source_homological_degree, int source_degree);

FreeModuleHomomorphism *ResolutionHomomorphism_getMap(ResolutionHomomorphism *f, uint homological_degree);


typedef struct {
    Resolution *resolution;
    Resolution *unit_resolution;
    uint max_product_homological_degree;
    ResolutionHomomorphism ****chain_maps_to_trivial_module_resolution;
} ResolutionWithMapsToUnitResolution;


ResolutionWithMapsToUnitResolution *ResolutionWithMapsToUnitResolution_construct(Resolution *res, Resolution *unit_res, uint max_homological_degree);
void ResolutionWithMapsToUnitResolution_extendMaps(ResolutionWithMapsToUnitResolution *res, uint homological_degree, int internal_degree);

#endif // CSTEENROD_RESOLUTION_HOMOMORPHISM_H