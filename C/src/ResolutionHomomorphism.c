#include <assert.h>
#include <limits.h>
#include <stdio.h>

#include "ResolutionHomomorphism.h"

ResolutionHomomorphism *ResolutionHomomorphism_construct(
    Resolution *source, Resolution *target, 
    uint homological_degree_shift, int internal_degree_shift
){
    assert(homological_degree_shift < target->max_homological_degree);
    uint max_hom_deg = source->max_homological_degree - homological_degree_shift;
    if(max_hom_deg > target->max_homological_degree){
        max_hom_deg = target->max_homological_degree;
    }
    int max_degree = source->max_internal_degree;
    if(target->max_internal_degree < max_degree){
        max_degree = target->max_internal_degree;
    }
    ResolutionHomomorphism *result = malloc(
        sizeof(ResolutionHomomorphism)
        + max_hom_deg * sizeof(FreeModuleHomomorphism*)
        + max_hom_deg * sizeof(int)
    );
    result->source = source;
    result->target = target;
    result->homological_degree_shift = homological_degree_shift;
    result->internal_degree_shift = internal_degree_shift;
    result->maps = (FreeModuleHomomorphism**)(result + 1);
    result->computed_degree = (int*)(result->maps + max_hom_deg);
    result->max_degree = max_degree;
    result->max_homological_degree = max_hom_deg;
    for(uint i = 0; i < max_hom_deg; i++){
        result->computed_degree[i] = result->source->min_internal_degree;   
    }
    result->maps[0] = FreeModuleHomomorphism_construct(
        source->modules[homological_degree_shift + 1], (Module*)target->modules[1], -internal_degree_shift, max_degree
    );
    for(int i = result->source->min_internal_degree; i < result->source->max_internal_degree; i++){
        uint num_gens = FreeModule_getNumberOfGensInDegree(result->source->modules[homological_degree_shift + 1], i);
        FreeModuleHomomorphism_AllocateSpaceForNewGenerators(result->maps[0], i, num_gens);
    }
    return result;
}

void ResolutionHomomorphism_setBaseMap(ResolutionHomomorphism *f, int gen_degree, int gen_index, Vector *output){
    assert(f->computed_degree[0] <= gen_degree);
    assert(gen_degree - f->internal_degree_shift >= f->target->min_internal_degree);
    // assert(gen_degree - f->internal_degree_shift >= f->target->min_degree );
    // assert(Resolution_cycleQ(f->target, 0, gen_degree - f->internal_degree_shift, output));
    FreeModuleHomomorphism_setOutput(f->maps[0], gen_degree, gen_index, output);
}

void ResolutionHomomorphism_baseMapReady(ResolutionHomomorphism *f, int degree){
    f->computed_degree[0] = degree + 1;
}

void ResolutionHomomorphism_extend_step(ResolutionHomomorphism *f, uint input_homological_degree, int input_internal_degree);
void ResolutionHomomorphism_extend(ResolutionHomomorphism *f, uint source_homological_degree, int source_degree){
    assert(f->computed_degree[0] >= source_degree);
    assert(source_homological_degree >= f->homological_degree_shift);
    for(uint i = 1; i < source_homological_degree - f->homological_degree_shift; i++){
        if(f->computed_degree[i] <= source_degree){
            if(f->computed_degree[i] == f->source->min_internal_degree){
                f->maps[i] = FreeModuleHomomorphism_construct(
                    f->source->modules[i + f->homological_degree_shift + 1], 
                    (Module*)f->target->modules[i + 1],
                    -f->internal_degree_shift,
                    f->max_degree
                );
            }
            for(int j = f->computed_degree[i]; j <= source_degree; j++){
                ResolutionHomomorphism_extend_step(f, i + f->homological_degree_shift, j);
            }
        }
    }
}

void ResolutionHomomorphism_extend_step(ResolutionHomomorphism *f, uint input_homological_degree, int input_internal_degree){
    uint p = f->source->algebra->p;
    assert(input_homological_degree >= f->homological_degree_shift);
    uint output_homological_degree = input_homological_degree - f->homological_degree_shift;
    f->computed_degree[output_homological_degree] = input_internal_degree + 1;
    int output_internal_degree = input_internal_degree - f->internal_degree_shift;        
    FreeModuleHomomorphism *d_source = f->source->differentials[input_homological_degree + 1];
    FreeModuleHomomorphism *d_target = f->target->differentials[output_homological_degree + 1];
    FreeModuleHomomorphism *f_cur = f->maps[output_homological_degree];
    FreeModuleHomomorphism *f_prev = f->maps[output_homological_degree - 1];
    assert(d_source->source == f_cur->source);
    assert(d_source->target == (Module*)f_prev->source);
    assert((Module*)d_target->source == f_cur->target);
    assert(d_target->target == f_prev->target);
    uint num_gens = FreeModule_getNumberOfGensInDegree(f_cur->source, input_internal_degree);
    FreeModuleHomomorphism_AllocateSpaceForNewGenerators(f_cur, input_internal_degree, num_gens);
    if(Module_getDimension((Module*)f->target->modules[output_homological_degree + 1], output_internal_degree) == 0){
        return; // nothing to do -- target is empty
    }
    Subspace *d_target_image = f->target->differentials[output_homological_degree]->kernel[output_internal_degree];
    Matrix *d_quasi_inverse = d_target->coimage_to_image_isomorphism[output_internal_degree];
    assert(d_quasi_inverse != NULL);
    uint dx_dimension = Module_getDimension((Module*)f_prev->source, input_internal_degree);
    uint fdx_dimension = Module_getDimension((Module*)f_prev->target, output_internal_degree);
    uint fx_dimension = Module_getDimension((Module*)f_cur->target, output_internal_degree);
    // assert(fdx_dimension == d_target_image->kernel->columns);
    char memory[
        Vector_getSize(p, dx_dimension, 0) 
        + Vector_getSize(p, fdx_dimension, 0)
        + Vector_getSize(p, fx_dimension, 0)
    ];
    char *memory_ptr = memory;
    Vector *dx_vector = Vector_initialize(p, &memory_ptr, dx_dimension, 0);
    Vector *fdx_vector = Vector_initialize(p, &memory_ptr, fdx_dimension, 0);
    Vector *fx_vector = Vector_initialize(p, &memory_ptr, fx_dimension, 0);
    for(uint k = 0; k < num_gens; k++){
        FreeModuleHomomorphism_applyToGenerator(d_source, dx_vector, 1, input_internal_degree, k);
        FreeModuleHomomorphism_apply(f_prev, fdx_vector, 1, input_internal_degree, dx_vector);
        Matrix_quasiInverse_apply(fx_vector, d_target_image, d_quasi_inverse, fdx_vector);
        FreeModuleHomomorphism_setOutput(f_cur, input_internal_degree, k, fx_vector);
        Vector_setToZero(dx_vector);
        Vector_setToZero(fdx_vector);
        Vector_setToZero(fx_vector);
    }
}

FreeModuleHomomorphism *ResolutionHomomorphism_getMap(ResolutionHomomorphism *f, uint homological_degree){
    return f->maps[homological_degree - f->homological_degree_shift];
}

ResolutionWithMapsToUnitResolution *ResolutionWithMapsToUnitResolution_construct(Resolution *res, Resolution *unit_res, uint max_homological_degree){
    size_t size = sizeof(ResolutionWithMapsToUnitResolution);
    size += sizeof(ResolutionHomomorphism***) * res->max_homological_degree;
    size += sizeof(ResolutionHomomorphism**) * (res->max_internal_degree - res->min_internal_degree);
    char *memory = malloc(size);
    ResolutionWithMapsToUnitResolution *result = (ResolutionWithMapsToUnitResolution *)memory;
    memory += sizeof(ResolutionWithMapsToUnitResolution);
    result->resolution = res;
    result->max_product_homological_degree = max_homological_degree;
    result->unit_resolution = unit_res;
    result->chain_maps_to_trivial_module_resolution = (ResolutionHomomorphism****)memory;
    memory += sizeof(ResolutionHomomorphism***) * res->max_homological_degree;
    for(uint i=0; i<res->max_homological_degree; i++){
        result->chain_maps_to_trivial_module_resolution[i] = (ResolutionHomomorphism***)memory;
        memory += sizeof(ResolutionHomomorphism**);
    }
    return result;
}

void ResolutionWithMapsToUnitResolution_extendMaps(ResolutionWithMapsToUnitResolution *res, uint homological_degree, int internal_degree){
    uint max_hom_deg = homological_degree;
    uint shifted_degree = internal_degree - res->resolution->min_internal_degree;
    if(res->max_product_homological_degree < max_hom_deg){
        max_hom_deg = res->max_product_homological_degree;
    }
    for(uint i = 0; i < max_hom_deg; i++){
        uint hom_deg = homological_degree - i;
        uint num_gens = res->resolution->modules[hom_deg + 1]->number_of_generators_in_degree[shifted_degree];
        for(uint j = 0; j < num_gens; j++){
            ResolutionHomomorphism_extend_step(
                res->chain_maps_to_trivial_module_resolution[hom_deg][shifted_degree][j],
                i, internal_degree
            );
        }
    }
}