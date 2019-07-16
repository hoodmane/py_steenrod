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
    uint homological_degree;
    int internal_degree;
    uint index;
} Cocycle;

typedef struct {
    uint length;
    uint capacity;
    Cocycle *list;
} Cocycle_list;

typedef struct {
    Resolution *resolution;
    Resolution *unit_resolution;
    uint max_product_homological_degree;
    Cocycle_list product_list;
    ResolutionHomomorphism ****chain_maps_to_trivial_module_resolution;
} ResolutionWithChainMaps;


ResolutionWithChainMaps *ResolutionWithChainMaps_construct(Resolution *res, Resolution *unit_res, uint number_of_products);
void ResolutionWithChainMaps_addProduct(ResolutionWithChainMaps *res_with_maps, uint homological_degree, int degree, uint index);


void ResolutionWithChainMaps_extendMaps(ResolutionWithChainMaps *res_with_maps, uint homological_degree, int internal_degree);
void ResolutionWithChainMaps_computeProduct(
    ResolutionWithChainMaps *res_with_maps, 
    uint elt_hom_deg, int elt_deg, uint elt_idx,
    uint source_hom_deg, int source_deg, uint source_idx
);

#endif // CSTEENROD_RESOLUTION_HOMOMORPHISM_H