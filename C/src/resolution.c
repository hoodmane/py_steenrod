//
// Created by Hood on 5/20/2019.
//

// TODO: split up this file a bit.

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "parson.h"
#include "Resolution.h"


#define max(a,b) ((a) > (b) ? (a) : (b))
void Resolution_generateOldKernelAndComputeNewKernel(Resolution *resolution, uint homological_degree, int degree);

void printCallback(Resolution * res, uint homological_degree, int degree, uint num_gens);

void addClass_doNothing(uint hom_deg __attribute__((unused)), int int_deg __attribute__((unused)), char *cocycle_name __attribute__((unused))){
    // printf("hom_deg: %d, int_deg: %d\n", hom_deg, int_deg);
}

void addStructline_doNothing(
    char *type __attribute__((unused)),
    uint source_hom_deg __attribute__((unused)), int source_int_deg __attribute__((unused)), uint source_idx __attribute__((unused)), 
    uint target_hom_deg __attribute__((unused)), int target_int_deg __attribute__((unused)), uint target_idx __attribute__((unused))
){
    // printf("\n\nhom_deg: %d, mod_deg: %d, num_gens: %d\n", homological_degree, degree, num_gens);
}

// Getter for javascript
FreeModuleHomomorphism *Resolution_getDifferential(Resolution *resolution, uint homological_degree){
    return resolution->differentials[homological_degree + 1];   
}

Resolution * Resolution_construct(
    FiniteDimensionalModule *module, 
    uint max_homological_degree,
    void (*addClass)(uint hom_deg, int int_deg, char *cocycle_name),
    void (*addStructline)(
        char *type,
        uint source_hom_deg, int source_int_deg, uint source_idx, 
        uint target_hom_deg, int target_int_deg, uint target_idx
    )    
){
    printf("Module name: %s\n", module->module.name);
    int min_degree = module->module.min_degree;
    // The 0th index in "modules" and "differentials" is the module we're resolving.
    // We're constructing an augmented resolution. In math this first module would be placed in index -1.
    // Anyways, we are working with the range -1 to max_degree which has max degree - (-1) + 1 = max_degree + 2 entries.
    Resolution *res = malloc(
        sizeof(Resolution)
        + (max_homological_degree + 1) * sizeof(FreeModule*)             // modules  
        + (max_homological_degree + 1) * sizeof(FreeModuleHomomorphism*) // differentials
        + (max_homological_degree + 1) * sizeof(int) // internal_degree_to_resolution_stage
    );
    res->modules = (FreeModule**)(res + 1);
    res->differentials = 
        (FreeModuleHomomorphism**)(res->modules + (max_homological_degree + 1));
    res->internal_degree_to_resolution_stage = (int*)(res->differentials + max_homological_degree + 1);
    memset(res->internal_degree_to_resolution_stage, 0, (max_homological_degree + 1)*sizeof(int));
    
    res->module = (Module*)module;
    res->algebra= module->module.algebra;
    res->computed_homological_degree = 0;
    res->max_homological_degree = max_homological_degree;
    res->modules[0] = (FreeModule*)module;
    res->differentials[0] = FreeModuleHomomorphism_construct((FreeModule*)module, NULL, 0, max_homological_degree);
    if(addClass == NULL){
        addClass = addClass_doNothing;
    }
    if(addStructline == NULL){
        addStructline = addStructline_doNothing;
    }    
    res->addClass = addClass;
    res->addStructline = addStructline;
    return res;
}

void Resolution_free(Resolution *res){
    if(res == NULL){
        return;
    }
    FreeModuleHomomorphism_free(res->differentials[0]);
    for(uint i = 0; i < res->max_homological_degree; i++){
        FreeModuleHomomorphism_free(res->differentials[i + 1]);
        FreeModule_free(res->modules[i + 1]);
    }
    free(res);
}

void Resolution_step(Resolution *res, uint homological_degree, int degree){
    uint shifted_degree = degree - res->module->min_degree;
    if(homological_degree == 0){
        FreeModuleHomomorphism *dminus1 = res->differentials[0];
        uint module_dim = Module_getDimension(res->module, degree);
        dminus1->kernel[shifted_degree] = Subspace_construct(res->algebra->p, module_dim, module_dim);
        for(uint j = 0; j < module_dim; j++){
            dminus1->kernel[shifted_degree]->column_to_pivot_row[j] = j;
            Vector_setEntry(dminus1->kernel[shifted_degree]->matrix->vectors[j], j, 1);
        }
    }
    // Construct new FreeModule.
    if(homological_degree >= res->computed_homological_degree){
        uint max_internal_degree = res->module->min_degree + res->max_homological_degree;
        res->computed_homological_degree ++;
        res->modules[homological_degree + 1] = 
            FreeModule_construct(res->algebra, res->module->min_degree, max_internal_degree);
        res->differentials[homological_degree + 1] =
            FreeModuleHomomorphism_construct(
                res->modules[homological_degree + 1], 
                (Module*)res->modules[homological_degree],
                0,
                max_internal_degree
            );
    }
    // Do the work
    // printf("(%d, %d)\n", homological_degree, degree);
    Resolution_generateOldKernelAndComputeNewKernel(res, homological_degree, degree);    
}

void Resolution_computeFiltrationOneProducts(Resolution *res, uint homological_degree, int degree, uint source_idx){
    FreeModuleHomomorphism *d = res->differentials[homological_degree + 1];
    FreeModule *T = (FreeModule*)d->target;        
    Vector *dx = d->outputs[degree - res->module->min_degree][source_idx];
    FiltrationOneProduct_list product_list = res->algebra->product_list; 
    for(uint j = 0; j < product_list.length; j++){
        char *op_name = product_list.list[j].type;
        uint op_degree = product_list.list[j].degree;
        uint op_index = product_list.list[j].index;
        int gen_degree = degree - op_degree;

        if(gen_degree < res->module->min_degree){
            break;
        }

        uint num_target_generators = T->number_of_generators_in_degree[gen_degree - res->module->min_degree];
        for(uint target = 0; target < num_target_generators; target++){
            uint vector_idx = FreeModule_operationGeneratorToIndex(res->modules[homological_degree], op_degree, op_index, gen_degree, target);
            if(vector_idx >= dx->dimension){
                printf("Out of bounds index when computing product:\n");
                printf("  ==  degree: %d, hom_deg: %d, dim: %d, idx: %d\n", degree, homological_degree, dx->dimension, vector_idx);
            } else {
                // printf("hom_deg: %d, deg: %d, source_idx: %d, op_deg: %d, entry: %d\n", homological_degree, degree, source_idx, op_degree, Vector_getEntry(dx, vector_idx));
                if(Vector_getEntry(dx, vector_idx) != 0){
                    // There was a product!
                    res->addStructline(op_name, homological_degree - 1, gen_degree, target, homological_degree, degree, source_idx);
                }
            }
        }
    }
}

//   It's assumed that before this function runs, we should have a complex whose homology H^{i,j} is 0 when 
//   i <= homological_degree, j <= degree and at least one of these inequalities is strict.
//   The output is a complex that's also exact in degree i, j.
//   Preconditions:
//      resolution->differentials[homological_degree]->kernel should contain the kernel of the previous differential
//      if homological_degree == 0, the kernel should be everything in the module.
void Resolution_generateOldKernelAndComputeNewKernel(Resolution *resolution, uint homological_degree, int degree){
    assert(degree >= (int)homological_degree + resolution->module->min_degree);
    uint shifted_degree = degree - resolution->module->min_degree;
    assert(shifted_degree < resolution->max_homological_degree);
    uint p = resolution->algebra->p;
    resolution->internal_degree_to_resolution_stage[shifted_degree] ++;
    FreeModuleHomomorphism *current_differential  = resolution->differentials[homological_degree + 1];
    uint source_dimension = Module_getDimension((Module*)current_differential->source, degree);
    uint target_dimension = Module_getDimension(current_differential->target, degree);
    assert(current_differential->source->number_of_generators_in_degree[shifted_degree] == 0);
    // The Homomorphism matrix has size source_dimension x target_dimension, but we are going to augment it with an
    // identity matrix so that gives a matrix with dimensions source_dimension x (target_dimension + source_dimension).
    // Later we're going to write into this same matrix an isomorphism source/image + new vectors --> kernel
    // This has size target_dimension x (2*target_dimension).
    // This latter matrix may be used to find a preimage of an element under the differential.

    // Pad the target dimension so that it ends in an aligned position.
    uint first_source_index = Vector_getPaddedDimension(p, target_dimension, 0);    
    uint rows = max(source_dimension, target_dimension);
    uint columns = first_source_index + source_dimension + rows;
    char matrix_memory[Matrix_getSize(p, rows, columns)];
    Matrix *matrix = Matrix_initialize(matrix_memory, p, rows, columns);
    
    int pivots[matrix->columns];
    FreeModuleHomomorphism_computeKernel(current_differential, matrix, pivots, degree);
    Subspace *kernel = current_differential->kernel[shifted_degree];

    // Now add generators to hit kernel of previous differential. 
    Subspace *previous_cycles = resolution->differentials[homological_degree]->kernel[shifted_degree];
    uint first_new_row = source_dimension - kernel->matrix->rows;
    // We stored the kernel rows somewhere else so we're going to write over them.
    // Add new free module generators to hit basis for previous kernel
    uint new_generators = Matrix_extendImage(matrix, first_source_index, source_dimension, first_new_row, pivots, previous_cycles);
    FreeModule_addGenerators(current_differential->source, degree, new_generators);
    FreeModuleHomomorphism_addGeneratorsFromMatrixRows(current_differential, degree, matrix, first_new_row, new_generators);

    // The part of the matrix that contains interesting information is occupied_rows x (target_dimension + source_dimension + kernel_size).
    // Allocate a matrix coimage_to_image with these dimensions.
    uint coimage_to_image_rows = first_new_row + new_generators;
    source_dimension += new_generators;
    Matrix *coimage_to_image = Matrix_construct(p, coimage_to_image_rows, source_dimension);
    current_differential->coimage_to_image_isomorphism[shifted_degree] = coimage_to_image;

    int useless_pivots[matrix->columns];
    rowReduce(matrix, useless_pivots, 0, 0);

    // Copy matrix contents to coimage_to_image
    char slice_memory[Vector_getSize(p, 0, 0)]; 
    char *slice_ptr = slice_memory;
    Vector *slice = Vector_initialize(p, &slice_ptr, 0, 0);
    for(uint i = 0; i < coimage_to_image_rows; i++) {
        Vector_slice(slice, matrix->vectors[i], first_source_index, first_source_index + source_dimension);
        Vector_assign(coimage_to_image->vectors[i], slice);
    }
}

uint Resolution_gradedDimensionString(char *buffer, Resolution *resolution){
    uint len = 0;
    len += sprintf(buffer + len, "[\n");
    for(uint i = 0; i < resolution->max_homological_degree; i++){
        len += sprintf(buffer + len, "  ");
        len += array_toString(buffer + len, 
                    &resolution->modules[i+1]->number_of_generators_in_degree[i], 
                    resolution->max_homological_degree - i
                );
        len += sprintf(buffer + len, ",\n");
    }
    len -= 2;
    len += sprintf(buffer + len, "\n]\n");
    return len;
}

uint Resolution_numberOfGensInDegree(Resolution *res, uint homological_degree, int internal_degree){
    return FreeModule_getNumberOfGensInDegree(res->modules[homological_degree + 1], internal_degree);
}

bool Resolution_cycleQ(Resolution *res, uint homological_degree, int internal_degree, Vector *element){
    if(homological_degree == 0){
        return true;
    }
    char buffer[2000];
    printf("hom_deg: %d, int_deg: %d\n", homological_degree, internal_degree);
    FreeModule_element_toJSONString(buffer, res->modules[homological_degree+1], internal_degree, element);
    printf("Element: %s\n", buffer);
    uint p = res->algebra->p;
    uint output_dimension = Module_getDimension((Module*)res->modules[homological_degree], internal_degree);
    char memory[Vector_getSize(p, output_dimension, 0)];
    char *memory_ptr = memory;
    Vector *output = Vector_initialize(p, &memory_ptr, output_dimension, 0);
    FreeModuleHomomorphism_apply(res->differentials[homological_degree+1], output, 1, internal_degree, element);
    FreeModule_element_toJSONString(buffer, res->modules[homological_degree], internal_degree, output);
    printf("Output: %s\n", buffer);
    return Vector_zeroQ(output);
}


SerializedResolution *Resolution_serialize(Resolution *res){
    JSON_Value *root_value = json_value_init_object();
    JSON_Object *root_object = json_object(root_value);
    json_object_set_number(root_object, "p", res->algebra->p);
    json_object_set_string(root_object, "algebra", res->algebra->name);
    json_object_set_string(root_object, "module", res->module->name);
    json_object_set_number(root_object, "max_homological_degree", res->computed_homological_degree);


    // Max computed degree
    JSON_Value *max_computed_degree_value = json_value_init_array();
    JSON_Array *max_computed_degree_array = json_array(max_computed_degree_value);
    for(uint i = 0; i < res->computed_homological_degree; i++){
        json_array_append_number(max_computed_degree_array, res->differentials[i+1]->max_computed_degree);
    }
    json_object_set_value(root_object, "max_computed_degrees", max_computed_degree_value);

    // Number of generators of each free module
    JSON_Value *dimensions_value = json_value_init_array();
    JSON_Array *dimensions_array = json_array(dimensions_value);
    for(uint i = 0; i < res->computed_homological_degree; i++){
        FreeModule *M = res->modules[i+1];
        FreeModuleHomomorphism *f = res->differentials[i+1];
        JSON_Value *current_module_dimensions_value = json_value_init_array();
        JSON_Array *current_module_dimensions_array = json_array(current_module_dimensions_value);
        for(int j = 0; j < f->max_computed_degree - M->module.min_degree; j++){
            json_array_append_number(current_module_dimensions_array, M->number_of_generators_in_degree[j]);
        }
        json_array_append_value(dimensions_array, current_module_dimensions_value);
    }
    json_object_set_value(root_object, "dimensions", dimensions_value);


    
    // differential "outputs" field
    JSON_Value *differentials_value = json_value_init_array();
    JSON_Array *differentials_array = json_array(differentials_value);
    size_t binary_size = 0;
    for(uint i=0; i<res->computed_homological_degree; i++){
        JSON_Value *ds1_value = json_value_init_array();
        JSON_Array *ds1_array = json_array(ds1_value);
        FreeModuleHomomorphism *f = res->differentials[i+1];
        for(int j=0; j < f->max_computed_degree; j++){
            if(f->coimage_to_image_isomorphism[j] == NULL){
                json_array_append_null(ds1_array);
            } else {
                Subspace *ker = f->kernel[j];
                Matrix *M = f->coimage_to_image_isomorphism[j];
                JSON_Value *matrix_value = json_value_init_object();
                JSON_Object *matrix_object = json_object(matrix_value);
                
                size_t vectors_size = 0;
                for(uint k=0; k<f->source->number_of_generators_in_degree[j]; k++){
                    vectors_size += f->outputs[j][k]->size;
                }
                json_object_set_number(matrix_object, "output_vectors_address", binary_size);                
                json_object_set_number(matrix_object, "output_vectors_size", vectors_size);
                binary_size += vectors_size;

                size_t kernel_size = Subspace_getSize(ker->matrix->p, ker->matrix->rows, ker->matrix->columns);
                json_object_set_number(matrix_object, "kernel_rows", ker->matrix->rows);
                json_object_set_number(matrix_object, "kernel_columns", ker->matrix->columns);
                json_object_set_number(matrix_object, "kernel_size", kernel_size);
                json_object_set_number(matrix_object, "kernel_address", binary_size); 
                binary_size += kernel_size;
                
                size_t matrix_size = Matrix_getSize(M->p, M->rows, M->columns);
                json_object_set_number(matrix_object, "preimage_rows", M->rows);
                json_object_set_number(matrix_object, "preimage_columns", M->rows);
                json_object_set_number(matrix_object, "preimage_size", matrix_size);
                json_object_set_number(matrix_object, "preimage_address", binary_size);
                binary_size += matrix_size;

                json_array_append_value(ds1_array, matrix_value);
            }            
        }
        json_array_append_value(differentials_array, ds1_value);
    }
    json_object_set_value(root_object, "differential_data", differentials_value);
    SerializedResolution *result = malloc(sizeof(SerializedResolution) + binary_size);
    result->json_data = json_serialize_to_string(root_value);
    result->json_size = strlen(result->json_data);
    char *binary_data = (char *)(result + 1);
    result->binary_size = binary_size;   
    result->binary_data = binary_data;

    for(uint i=0; i<res->computed_homological_degree; i++){
        FreeModuleHomomorphism *f = res->differentials[i+1];
        for(int j=0; j < f->max_computed_degree; j++){
            if(f->coimage_to_image_isomorphism[j] == NULL){
                continue;
            }            
            for(uint k=0; k<f->source->number_of_generators_in_degree[j]; k++){
                  Vector_serialize(&binary_data, f->outputs[j][k]);
            }
            Subspace_serialize(&binary_data, f->kernel[j]);
            Matrix_serialize(&binary_data, f->coimage_to_image_isomorphism[j]);
        }
    }
    return result;
}

Resolution *Resolution_deserialize(
    FiniteDimensionalModule *module, SerializedResolution *sres,
    void (*addClass)(uint hom_deg, int int_deg, char *cocycle_name),
    void (*addStructline)(
        char *type,
        uint source_hom_deg, int source_int_deg, uint source_idx, 
        uint target_hom_deg, int target_int_deg, uint target_idx
    )    
){
    JSON_Value *root_value = json_parse_string(sres->json_data);
    JSON_Object *root_object = json_object(root_value);
    uint max_homological_degree = json_object_get_number(root_object, "max_homological_degree");
    JSON_Array *dimensions = json_object_get_array(root_object, "dimensions");
    JSON_Array *max_computed_degrees = json_object_get_array(root_object, "max_computed_degrees");
    Resolution *res = Resolution_construct(module, max_homological_degree, addClass, addStructline);
    // Construct new FreeModule.
    for(uint homological_degree = 0; homological_degree < max_homological_degree; homological_degree++){
        res->computed_homological_degree ++;
        FreeModule *F = 
            FreeModule_construct(res->algebra, res->module->min_degree, res->module->min_degree + res->max_homological_degree);
        res->modules[homological_degree + 1] = F;
        JSON_Array *cur_dimensions = json_array_get_array(dimensions, homological_degree);
        for(uint i = 0; i < json_array_get_count(cur_dimensions); i++){
            F->number_of_generators_in_degree[i] = json_array_get_number(cur_dimensions, i);
            FreeModule_ConstructBlockOffsetTable(F, i);
        }
        FreeModuleHomomorphism *f =
            FreeModuleHomomorphism_construct(
                res->modules[homological_degree + 1], 
                (Module*)res->modules[homological_degree],
                0,
                res->max_homological_degree + res->module->min_degree
            );
        res->differentials[homological_degree + 1] = f;
        int max_computed_degree = json_array_get_number(max_computed_degrees, homological_degree); 
        for(int i=0; i < max_computed_degree; i++){
            FreeModuleHomomorphism_AllocateSpaceForNewGenerators(f, i, f->source->number_of_generators_in_degree[i]);
        }
    }
    // for(uint i = 0; i < homology_dimension; i++){
    //     FreeModuleHomomorphism_setOutput(current_differential, degree, i, slice);
    // }
    JSON_Array *matrices_outer_array = json_object_get_array(root_object, "differential_data");
    char *binary_data = sres->binary_data;
    for(uint i=0; i<res->computed_homological_degree; i++){
        JSON_Array *matrices_inner_array = json_array_get_array(matrices_outer_array, i);        
        FreeModuleHomomorphism *f = res->differentials[i+1];
        for(int j=0; j < f->max_computed_degree; j++){
            JSON_Object *obj = json_array_get_object(matrices_inner_array, j);
            if(obj == NULL){
                continue;
            }
            for(uint k=0; k<f->source->number_of_generators_in_degree[j]; k++){
                f->outputs[j][k] = Vector_deserialize(res->algebra->p, &binary_data);
            }
            f->kernel[j] = Subspace_deserialize(&binary_data);
            f->coimage_to_image_isomorphism[j] = Matrix_deserialize(&binary_data);
        }
    }
    return res;
}

size_t SerializedResolution_getJSONSize(SerializedResolution *sres){
    return sres->json_size;
}

char *SerializedResolution_getJSONData(SerializedResolution *sres){
    return sres->json_data;
}

size_t SerializedResolution_getBinarySize(SerializedResolution *sres){
    return sres->binary_size;
}

char *SerializedResolution_getBinaryData(SerializedResolution *sres){
    return sres->binary_data;
}