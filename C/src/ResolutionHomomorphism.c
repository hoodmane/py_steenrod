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
    int max_degree = source->max_homological_degree + source->module->min_degree;
    int target_max_degree = target->max_homological_degree + target->module->min_degree;
    if(target_max_degree < max_degree){
        max_degree = target_max_degree;
    }
    size_t size = sizeof(ResolutionHomomorphism)
                + max_hom_deg * sizeof(FreeModuleHomomorphism*)
                + max_hom_deg * sizeof(int);
    ResolutionHomomorphism *result = malloc(size);
    result->source = source;
    result->target = target;
    result->homological_degree_shift = homological_degree_shift;
    result->internal_degree_shift = internal_degree_shift;
    result->maps = (FreeModuleHomomorphism**)(result + 1);
    result->computed_degree = (int*)(result->maps + max_hom_deg);
    assert((char*)(result->computed_degree + max_hom_deg) == (char*)result + size); // Did we use the right amount of memory?
    result->max_degree = max_degree;
    result->max_homological_degree = max_hom_deg;
    int min_internal_degree = result->source->module->min_degree;
    int max_internal_degree = result->source->module->min_degree + result->source->max_homological_degree;
    for(uint i = 0; i < max_hom_deg; i++){
        result->computed_degree[i] = min_internal_degree;   
    }
    result->maps[0] = FreeModuleHomomorphism_construct(
        source->modules[homological_degree_shift + 1], (Module*)target->modules[1], -internal_degree_shift, max_degree
    );
    ResolutionHomomorphism_extendBaseMap(result, result->maps[0]->source->module.max_computed_degree);
    return result;
}

void ResolutionHomomorphism_extendBaseMap(ResolutionHomomorphism *f, int degree){
    for(int j = f->maps[0]->max_computed_degree; j < degree; j++){
        uint num_gens = FreeModule_getNumberOfGensInDegree(f->maps[0]->source, j);
        FreeModuleHomomorphism_AllocateSpaceForNewGenerators(f->maps[0], j, num_gens);
    }
}

void ResolutionHomomorphism_setBaseMap(ResolutionHomomorphism *f, int gen_degree, int gen_index, Vector *output){
    assert(f->computed_degree[0] <= gen_degree);
    assert(gen_degree - f->internal_degree_shift >= f->target->module->min_degree);
    // assert(gen_degree - f->internal_degree_shift >= f->target->min_degree );
    // assert(Resolution_cycleQ(f->target, 0, gen_degree - f->internal_degree_shift, output));
    FreeModuleHomomorphism_setOutput(f->maps[0], gen_degree, gen_index, output);
}

void ResolutionHomomorphism_baseMapReady(ResolutionHomomorphism *f, int degree){
    f->computed_degree[0] = degree + 1;
}

void ResolutionHomomorphism_extend_step(ResolutionHomomorphism *f, uint input_homological_degree, int input_internal_degree);
void ResolutionHomomorphism_extend(ResolutionHomomorphism *f, uint source_homological_degree, int source_degree){
    // printf("      f: %llx, hd: %d, id: %d\n",(uint64)f,source_homological_degree, source_degree);
    assert(f->computed_degree[0] >= source_degree);
    if(source_homological_degree < f->homological_degree_shift){
        return;
    }
    ResolutionHomomorphism_extendBaseMap(f, source_degree);

    for(uint i = 1; i <= source_homological_degree - f->homological_degree_shift; i++){
        if(f->computed_degree[i] == f->source->module->min_degree){
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

ResolutionWithChainMaps *ResolutionWithChainMaps_construct(Resolution *res, Resolution *unit_res, uint number_of_products){
    size_t size = sizeof(ResolutionWithChainMaps);
    size += sizeof(Cocycle) * number_of_products;
    size += sizeof(ResolutionHomomorphism***) * res->max_homological_degree;
    size += sizeof(ResolutionHomomorphism**) * res->max_homological_degree * res->max_homological_degree;
    char *memory = malloc(size);
    ResolutionWithChainMaps *result = (ResolutionWithChainMaps *)memory;
    memory += sizeof(ResolutionWithChainMaps);
    result->resolution = res;
    result->max_product_homological_degree = 0;
    result->unit_resolution = unit_res;
    result->product_list.capacity = number_of_products;
    result->product_list.length = 0;
    result->product_list.list = (Cocycle*) memory;
    memory += number_of_products * sizeof(Cocycle);
    result->chain_maps_to_trivial_module_resolution = (ResolutionHomomorphism****)memory;
    memory += sizeof(ResolutionHomomorphism***) * res->max_homological_degree;
    for(uint i=0; i<res->max_homological_degree; i++){
        result->chain_maps_to_trivial_module_resolution[i] = (ResolutionHomomorphism***)memory;
        memory += sizeof(ResolutionHomomorphism**) * res->max_homological_degree;
    }
    assert((char*)result + size == memory);
    return result;
}

void ResolutionWithChainMaps_addProduct(ResolutionWithChainMaps *res_with_maps, uint homological_degree, int degree, uint index){
    assert(res_with_maps->product_list.length < res_with_maps->product_list.capacity);
    if(homological_degree > res_with_maps->max_product_homological_degree){
        res_with_maps->max_product_homological_degree = homological_degree;
    }
    Cocycle *c = &res_with_maps->product_list.list[res_with_maps->product_list.length];
    res_with_maps->product_list.length ++;
    c->homological_degree = homological_degree;
    c->internal_degree = degree;
    c->index = index;
}

void ResolutionWithChainMaps_extendMaps(ResolutionWithChainMaps *res_with_maps, uint homological_degree, int internal_degree){
    assert(res_with_maps->product_list.capacity == res_with_maps->product_list.length);
    if(res_with_maps->max_product_homological_degree == 0){
        return;
    }
    // printf("   ems: (%d, %d)\n", homological_degree,internal_degree);
    uint max_hom_deg = homological_degree;
    uint shifted_degree = internal_degree - res_with_maps->resolution->module->min_degree;
    if(res_with_maps->max_product_homological_degree < max_hom_deg){
        max_hom_deg = res_with_maps->max_product_homological_degree;
    }
    uint p = res_with_maps->resolution->algebra->p;
    char memory[Vector_getSize(p, 1, 0)];
    char *memory_ptr = memory;
    Vector *unit_vector = Vector_initialize(p, &memory_ptr, 1, 0);
    Vector_setEntry(unit_vector, 0, 1);
    uint hom_deg = homological_degree;
    uint num_gens = FreeModule_getNumberOfGensInDegree(res_with_maps->resolution->modules[hom_deg + 1], internal_degree);
    if(num_gens > 0){
        ResolutionHomomorphism*** table = &res_with_maps->chain_maps_to_trivial_module_resolution[hom_deg][shifted_degree];
        *table = malloc(num_gens * sizeof(ResolutionHomomorphism*));
        for(uint j=0; j<num_gens; j++){
            ResolutionHomomorphism *f = ResolutionHomomorphism_construct(res_with_maps->resolution, res_with_maps->unit_resolution, homological_degree, internal_degree);
            ResolutionHomomorphism_extendBaseMap(f, internal_degree + 1);
            ResolutionHomomorphism_setBaseMap(f, internal_degree, j, unit_vector);
            ResolutionHomomorphism_baseMapReady(f, f->maps[0]->source->module.max_computed_degree);
            (*table)[j] = f;
        }
    }

    int min_degree = res_with_maps->resolution->module->min_degree;
    for(uint i = 0; i <= max_hom_deg; i++){
        for(int j = min_degree; j <= internal_degree; j++){
            hom_deg = homological_degree - i;
            uint num_gens = FreeModule_getNumberOfGensInDegree(res_with_maps->resolution->modules[hom_deg + 1], j);
            for(uint k = 0; k < num_gens; k++){
                // printf("      cocyc (%d, %d, %d) to (%d, %d) \n", hom_deg, j, k,  i, internal_degree);
                ResolutionHomomorphism *f = res_with_maps->chain_maps_to_trivial_module_resolution[hom_deg][j - min_degree][k];
                ResolutionHomomorphism_extendBaseMap(f, internal_degree + 1);
                f->computed_degree[0] = internal_degree + 1;                
                ResolutionHomomorphism_extend(f, homological_degree, internal_degree);
            }
        }
    }
}

void ResolutionWithChainMaps_computeProduct(
    ResolutionWithChainMaps *res_with_maps, 
    uint elt_hom_deg, int elt_deg, uint elt_idx,
    uint source_hom_deg, int source_deg, uint source_idx
){
    Resolution *res = res_with_maps->resolution;
    ResolutionHomomorphism *f = res_with_maps->chain_maps_to_trivial_module_resolution[source_hom_deg][source_deg][source_idx];
    uint target_hom_deg = source_hom_deg + elt_hom_deg;
    uint target_deg = source_deg + elt_deg;
    FreeModule *output_module = res->modules[elt_hom_deg + 1];
    uint output_gens = FreeModule_getNumberOfGensInDegree(output_module, elt_deg);
    Vector *output = Vector_construct(res->algebra->p, Module_getDimension((Module*)output_module,  elt_deg), 0);
    for(uint l = 0; l < Resolution_numberOfGensInDegree(res, target_hom_deg, target_deg); l++){
        FreeModuleHomomorphism_applyToGenerator(f->maps[elt_hom_deg], output, 1, target_deg, l);
        if(Vector_getEntry(output, output->dimension - output_gens + elt_idx) != 0){
            res->addStructline(
                source_hom_deg, source_deg, source_idx, 
                target_hom_deg, target_deg, l
            );
        }
    }    
}