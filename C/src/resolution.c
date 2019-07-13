//
// Created by Hood on 5/20/2019.
//

// TODO: split up this file a bit.

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "parson.h"
#include "combinatorics.h"
#include "Algebra.h"
#include "AdemAlgebra.h"
#include "MilnorAlgebra.h"
#include "Module.h"
#include "Resolution.h"
#include "ResolutionHomomorphism.h"

#define max(a,b) ((a) > (b) ? (a) : (b))
void Resolution_generateOldKernelAndComputeNewKernel(Resolution *resolution, uint homological_degree, int degree);

void printCallback(Resolution * res, uint homological_degree, int degree, uint num_gens);

void addClass_doNothing(uint hom_deg __attribute__((unused)), int int_deg __attribute__((unused)), char *cocycle_name __attribute__((unused))){
    // printf("hom_deg: %d, int_deg: %d\n", hom_deg, int_deg);
}

void addStructline_doNothing(
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
    int max_degree,
    void (*addClass)(uint hom_deg, int int_deg, char *cocycle_name),
    void (*addStructline)(
        uint source_hom_deg, int source_int_deg, uint source_idx, 
        uint target_hom_deg, int target_int_deg, uint target_idx
    )    
){
    printf("Module name: %s\n", module->module.name);
    int min_degree = module->module.min_degree;
    uint num_degrees = max_degree - min_degree;
    // The 0th index in "modules" and "differentials" is the module we're resolving.
    // We're constructing an augmented resolution. In math this first module would be placed in index -1.
    // Anyways, we are working with the range -1 to max_degree which has max degree - (-1) + 1 = max_degree + 2 entries.
    Resolution *res = malloc(
        sizeof(Resolution)
        + (num_degrees + 1) * sizeof(FreeModule*)
        + (num_degrees + 1) * sizeof(FreeModuleHomomorphism*) // 
        + (num_degrees + 1) * sizeof(int) // internal_degree_to_resolution_stage
    );
    res->modules = (FreeModule**)(res + 1);
    res->differentials = 
        (FreeModuleHomomorphism**)(res->modules + (num_degrees + 1));
    res->internal_degree_to_resolution_stage = (int*)(res->differentials + num_degrees + 1);
    memset(res->internal_degree_to_resolution_stage, 0, (num_degrees + 1)*sizeof(int));
    
    res->module = (Module*)module;
    res->algebra= module->module.algebra;
    res->max_homological_degree = 0;
    res->min_degree = res->module->min_degree;
    res->max_degree = max_degree;
    res->modules[0] = (FreeModule*)module;
    res->differentials[0] = FreeModuleHomomorphism_construct((FreeModule*)module, NULL, 0, max_degree);
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
    for(int i = 0; i < res->max_degree; i++){
        FreeModuleHomomorphism_free(res->differentials[i + 1]);
        FreeModule_free(res->modules[i + 1]);
    }
    free(res);
}

void Resolution_resolveThroughDegree(Resolution *res, int degree){
//     for(uint hom_deg = 0; hom_deg < degree; hom_deg++){
//        for(uint int_deg = hom_deg; int_deg <= degree; int_deg ++){
//             Resolution_step(res, hom_deg, int_deg);
//             printf("(%d, %d)\n", hom_deg, int_deg);
//        }
//    }
    for(int int_deg = res->min_degree; int_deg < degree; int_deg ++){
        for(int hom_deg = 0; hom_deg <= int_deg - res->min_degree; hom_deg++){           
            // printf("(%d, %d)\n", hom_deg, int_deg);
            Resolution_step(res, hom_deg, int_deg);
        }
    }
    // printf
    // Resolution_serialize(res);
    // for(int i = degree - 1 - res->min_degree; i >= 0; i--){
    //     printf("stage %*d: ", 2, i);
    //     array_print("%s\n", &res->modules[i+1]->number_of_generators_in_degree[i], degree - i - res->min_degree);
    // }   
}

void Resolution_computeFiltrationOneProducts(Resolution *res, uint homological_degree, int degree, uint source_idx);
void Resolution_step(Resolution *res, uint homological_degree, int degree){
    // printf("degree: %d, homological_degree: %d\n", degree, homological_degree);
    // Construct kernel -- say that it's everything.
    // We put the module itself in degree zero and we'll want to hit the whole thing.
    uint shifted_degree = degree - res->min_degree;
    if(homological_degree == 0){
        FreeModuleHomomorphism *dminus1 = res->differentials[0];
        uint module_dim = Module_getDimension(res->module, degree);
        dminus1->kernel[shifted_degree] = Kernel_construct(res->algebra->p, module_dim, module_dim);
        for(uint j = 0; j < module_dim; j++){
            dminus1->kernel[shifted_degree]->column_to_pivot_row[j] = j;
            Vector_setEntry(dminus1->kernel[shifted_degree]->kernel->matrix[j], j, 1);
        }
    }
    // Construct new FreeModule.
    if(homological_degree >= res->max_homological_degree){
        res->max_homological_degree ++;
        res->modules[homological_degree + 1] = 
            FreeModule_construct(res->algebra, res->min_degree, res->max_degree);
        res->differentials[homological_degree + 1] =
            FreeModuleHomomorphism_construct(
                res->modules[homological_degree + 1], 
                (Module*)res->modules[homological_degree],
                0,
                res->max_degree
            );
    }
    // Do the work
    // printf("(%d, %d)\n", homological_degree, degree);
    Resolution_generateOldKernelAndComputeNewKernel(res, homological_degree, degree);

    // Report the answers.
    // Classes:
    uint num_gens = res->modules[homological_degree + 1]->number_of_generators_in_degree[degree-res->min_degree];
    for(uint i=0; i < num_gens; i++){
        // printf("addClass(%d, %d)\n", homological_degree, degree);
        res->addClass(homological_degree, degree, "");
        if(homological_degree > 0){
            Resolution_computeFiltrationOneProducts(res, homological_degree, degree, i);
        }
    } 
    // Products. TODO: handle case distinction by primes.

}

void Resolution_computeFiltrationOneProducts(Resolution *res, uint homological_degree, int degree, uint source_idx){
    FreeModuleHomomorphism *d = res->differentials[homological_degree + 1];
    FreeModule *T = (FreeModule*)d->target;        
    Vector *dx = d->outputs[degree - res->min_degree][source_idx];
    FiltrationOneProductList *product_list = res->algebra->product_list; 
    for(uint j = 0; j < product_list->length; j++){
        uint op_degree = product_list->degrees[j];
        uint op_index = product_list->indices[j];
        int gen_degree = degree - op_degree;

        if(gen_degree < res->min_degree){
            break;
        }

        uint num_target_generators = T->number_of_generators_in_degree[gen_degree - res->min_degree];
        for(uint target = 0; target < num_target_generators; target++){
            uint vector_idx = FreeModule_operationGeneratorToIndex(res->modules[homological_degree], op_degree, op_index, gen_degree, target);
            if(vector_idx >= dx->dimension){
                printf("Out of bounds index when computing product:\n");
                printf("  ==  degree: %d, hom_deg: %d, dim: %d, idx: %d\n", degree, homological_degree, dx->dimension, vector_idx);
            } else {
                // printf("hom_deg: %d, deg: %d, source_idx: %d, op_deg: %d, entry: %d\n", homological_degree, degree, source_idx, op_degree, Vector_getEntry(dx, vector_idx));
                if(Vector_getEntry(dx, vector_idx) != 0){
                    // There was a product!
                    res->addStructline(homological_degree - 1, gen_degree, target, homological_degree, degree, source_idx);
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
    // printf("degree: %d, homological_degree: %d, resolution->min_degree: %d\n",degree, homological_degree, resolution->min_degree);
    // printf("(%d, %d)\n", homological_degree, degree);
    assert(degree >= (int)homological_degree + resolution->min_degree);
    assert(degree < resolution->max_degree);
    uint shifted_degree = degree - resolution->module->min_degree;
    uint p = resolution->algebra->p;
    FreeModuleHomomorphism *current_differential  = resolution->differentials[homological_degree + 1];
    FreeModuleHomomorphism *previous_differential = resolution->differentials[homological_degree];
    FreeModule *source = current_differential->source;

    uint source_dimension = Module_getDimension(&current_differential->source->module, degree);
    uint target_dimension = Module_getDimension(current_differential->target, degree);
    assert(current_differential->source->number_of_generators_in_degree[shifted_degree] == 0);

    // The Homomorphism matrix has size source_dimension x target_dimension, but we are going to augment it with an
    // identity matrix so that gives a matrix with dimensions source_dimension x (target_dimension + source_dimension).
    // Later we're going to write into this same matrix an isomorphism source/image + new vectors --> kernel
    // This has size target_dimension x (2*target_dimension).
    // This latter matrix may be used to find a preimage of an element under the differential.
    // The information about the kernel is moved to a separate temporary kernel matrix and eventually can be found by looking at
    // the image of the next step in the resolution.

    // Pad the target dimension so that it ends in an aligned position.
    uint padded_target_dimension = Vector_getPaddedDimension(p, target_dimension, 0);
    uint rows = max(source_dimension, target_dimension);
    uint columns = padded_target_dimension + source_dimension + rows;
    char matrix_memory[Matrix_getSize(p, rows, columns)];
    Matrix *full_matrix = Matrix_initialize(matrix_memory, p, rows, columns);
    // For the first stage we just want the part of size source_dimension x (padded_target_dimension + source_dimension).
    // Slice matrix out of full_matrix.
    char slice_matrix_memory[Matrix_getSliceSize(p, source_dimension)];
    Matrix *matrix = Matrix_slice(full_matrix, slice_matrix_memory, 0, source_dimension, 0, padded_target_dimension + source_dimension);
    FreeModuleHomomorphism_getMatrix(current_differential, matrix, degree);

    // Write the identity matrix into the right block
    for(uint i = 0; i < source_dimension; i++){
        Vector_setEntry(matrix->matrix[i], padded_target_dimension + i, 1);
    }

    // Row reduce
    int column_to_pivot_row[matrix->columns];
    rowReduce(matrix, column_to_pivot_row, 0, 0);//target_dimension, padded_target_dimension);

    uint permutation[matrix->rows];
    Matrix_getRowPermutation(matrix, permutation);
    Matrix_applyRowPermutation(full_matrix, permutation, matrix->rows);

    // Stage 1: Find kernel of current differential
    // Locate first kernel row
    uint first_kernel_row = matrix->rows;
    for(uint i = padded_target_dimension; i < matrix->columns; i ++){
        if(column_to_pivot_row[i] != -1){
            first_kernel_row = column_to_pivot_row[i];
            break;
        }
    }
    // Every row after the first kernel row is also a kernel row, so now we know how big it is and can allocate space.
    uint kernel_dimension = matrix->rows - first_kernel_row;
    Kernel *kernel = Kernel_construct(p, kernel_dimension, source_dimension);
    // Write pivots into kernel
    for(uint i = 0; i < source_dimension; i++){
        // Turns -1 into some negative number... make sure to check <0 for no pivot in column...
        kernel->column_to_pivot_row[i] = column_to_pivot_row[i + padded_target_dimension] - first_kernel_row;
    }

    // Copy kernel matrix into kernel
    for(uint row = 0; row < kernel_dimension; row++){
        char slice_memory[Vector_getSize(p, 0, 0)];
        char *slice_ptr = slice_memory;
        Vector *slice = Vector_initialize(p, &slice_ptr, 0, 0);
        Vector_slice(slice, matrix->matrix[first_kernel_row + row], padded_target_dimension, padded_target_dimension + source_dimension);
        Vector_assign(kernel->kernel->matrix[row], slice);
    }
    current_differential->kernel[shifted_degree] = kernel;

    // Stage 2: Hit kernel of previous differential. 
    Kernel *previous_cycles = previous_differential->kernel[shifted_degree];
    // uint previous_cycle_dimension = previous_cycles->kernel->rows;
    // We no longer care about the kernel rows since we stored them somewhere else, 
    // so we're going to write over them.
    uint current_target_row = first_kernel_row;
    // Find basis of quotient previous_kernel/image and add new free module generators to hit
    uint homology_dimension = 0;
    for(uint i = 0; i < previous_cycles->kernel->columns; i++){
        if(column_to_pivot_row[i] < 0 && previous_cycles->column_to_pivot_row[i] >= 0){
            // Look up the cycle that we're missing and add a generator hitting it.
            int kernel_vector_row = previous_cycles->column_to_pivot_row[i];
            // assert(kernel_vector_row < previous_kernel->kernel->rows);
            Vector *new_image = previous_cycles->kernel->matrix[kernel_vector_row];
            // Stack allocate slice
            char slice_memory[Vector_getSize(p, 0, 0)];
            char *slice_ptr = slice_memory;
            Vector *slice = Vector_initialize(p, &slice_ptr, 0, 0);
            // Write new image to full_matrix
            Vector_slice(slice, full_matrix->matrix[current_target_row], 0, previous_cycles->kernel->columns);
            Vector_assign(slice, new_image);
            // Vector_print("    new_image: %s\n", slice);
            // Write elementary basis vector into right block of full_matrix
            Vector_slice(slice, full_matrix->matrix[current_target_row], padded_target_dimension, full_matrix->columns);
            Vector_setToZero(slice);
            Vector_setEntry(slice, source_dimension + homology_dimension, 1);
            current_target_row++;
            homology_dimension++;
        }
    }
    source->number_of_generators_in_degree[shifted_degree] = homology_dimension;
    // Copy the outputs, currently stored in the coimage_to_image matrix, to the FreeModuleHomomorphism outputs field.
    FreeModuleHomomorphism_AllocateSpaceForNewGenerators(current_differential, degree, homology_dimension);
    for(uint i = 0; i < homology_dimension; i++){
        char slice_memory[Vector_getSize(p, 0, 0)]; 
        char *slice_ptr = slice_memory;
        Vector *slice = Vector_initialize(p, &slice_ptr, 0, 0);
        Vector_slice(slice, full_matrix->matrix[first_kernel_row + i], 0, target_dimension);
        FreeModuleHomomorphism_setOutput(current_differential, degree, i, slice);
    }
    FreeModule_ConstructBlockOffsetTable(source, degree);

    // The part of the matrix that contains interesting information is current_target_row * (target_dimension + source_dimension + kernel_size).
    // Allocate a matrix coimage_to_image with these dimensions.
    uint coimage_to_image_rows = current_target_row;
    uint coimage_to_image_column_start = padded_target_dimension;
    uint coimage_to_image_columns = source_dimension + homology_dimension;
    Matrix *coimage_to_image = Matrix_construct(p, coimage_to_image_rows, coimage_to_image_columns);
    current_differential->coimage_to_image_isomorphism[shifted_degree] = coimage_to_image;

    int useless_pivot_row_info[full_matrix->columns];
    rowReduce(full_matrix, useless_pivot_row_info, 0, 0);

    // Copy matrix contents to coimage_to_image
    for(uint i = 0; i < coimage_to_image_rows; i++) {
        char slice_memory[Vector_getSize(p, 0, 0)]; 
        char *slice_ptr = slice_memory;
        Vector *slice = Vector_initialize(p, &slice_ptr, 0, 0);
        Vector_slice(slice, full_matrix->matrix[i], coimage_to_image_column_start, coimage_to_image_column_start + coimage_to_image_columns);
        Vector_assign(coimage_to_image->matrix[i], slice);
    }


    // TODO: assertion about contents of useless_pivot_row_info?
    // Should contain [0,1,2,3,...,n,-1,-1,-1,..., -1].
    resolution->internal_degree_to_resolution_stage[shifted_degree] ++;
}

uint Resolution_gradedDimensionString(char *buffer, Resolution *resolution){
    uint len = 0;
    len += sprintf(buffer + len, "[\n");
    for(int i = 0; i < resolution->max_degree; i++){
        len += sprintf(buffer + len, "  ");
        len += array_toString(buffer + len, 
                    &resolution->modules[i+1]->number_of_generators_in_degree[i], 
                    resolution->max_degree - i
                );
        len += sprintf(buffer + len, ",\n");
    }
    len -= 2;
    len += sprintf(buffer + len, "\n]\n");
    return len;
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
    json_object_set_number(root_object, "max_homological_degree", res->max_homological_degree);


    // Max computed degree
    JSON_Value *max_computed_degree_value = json_value_init_array();
    JSON_Array *max_computed_degree_array = json_array(max_computed_degree_value);
    for(uint i = 0; i < res->max_homological_degree; i++){
        json_array_append_number(max_computed_degree_array, res->differentials[i+1]->max_computed_degree);
    }
    json_object_set_value(root_object, "max_computed_degrees", max_computed_degree_value);

    // Number of generators of each free module
    JSON_Value *dimensions_value = json_value_init_array();
    JSON_Array *dimensions_array = json_array(dimensions_value);
    for(uint i = 0; i < res->max_homological_degree; i++){
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
    for(uint i=0; i<res->max_homological_degree; i++){
        JSON_Value *ds1_value = json_value_init_array();
        JSON_Array *ds1_array = json_array(ds1_value);
        FreeModuleHomomorphism *f = res->differentials[i+1];
        for(int j=0; j < f->max_computed_degree; j++){
            if(f->coimage_to_image_isomorphism[j] == NULL){
                json_array_append_null(ds1_array);
            } else {
                Kernel *ker = f->kernel[j];
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

                size_t kernel_size = Kernel_getSize(ker->kernel->p, ker->kernel->rows, ker->kernel->columns);
                json_object_set_number(matrix_object, "kernel_rows", ker->kernel->rows);
                json_object_set_number(matrix_object, "kernel_columns", ker->kernel->columns);
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

    for(uint i=0; i<res->max_homological_degree; i++){
        FreeModuleHomomorphism *f = res->differentials[i+1];
        for(int j=0; j < f->max_computed_degree; j++){
            if(f->coimage_to_image_isomorphism[j] == NULL){
                continue;
            }            
            for(uint k=0; k<f->source->number_of_generators_in_degree[j]; k++){
                  Vector_serialize(&binary_data, f->outputs[j][k]);
            }
            Kernel_serialize(&binary_data, f->kernel[j]);
            Matrix_serialize(&binary_data, f->coimage_to_image_isomorphism[j]);
        }
    }
    return result;
}

Resolution *Resolution_deserialize(
    FiniteDimensionalModule *module, SerializedResolution *sres,
    void (*addClass)(uint hom_deg, int int_deg, char *cocycle_name),
    void (*addStructline)(
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
        res->max_homological_degree ++;
        FreeModule *F = 
            FreeModule_construct(res->algebra, res->min_degree, res->max_degree);
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
                res->max_degree
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
    for(uint i=0; i<res->max_homological_degree; i++){
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
            f->kernel[j] = Kernel_deserialize(&binary_data);
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

/**/
int main(int argc, char *argv[]){
    // assert(false);
    uint p;
    uint degree;
    if(argc > 2){
        char *end;
        p = strtoul(argv[2], &end, 10);
    } else {
        p = 2;
    }

    if(argc > 3){
        char *end;
        degree = strtoul(argv[3], &end, 10);
    } else {
        degree = 20;
    }

    // uint degree = 70;
    bool generic = p!=2;
    initializePrime(p);

    // uint p_part[3] = {3,2,1};
    // uint p_part_length = 3;
    // uint p_part[2] = {2,1};
    // uint p_part_length = 2;
    // bool truncated = true;
    // Profile *P = NULL; // Profile_construct(generic, 0, NULL, p_part_length, p_part, truncated);
    // MilnorAlgebra *A = MilnorAlgebra_construct(p, generic, P);
    Algebra *algebra;
    if(argc > 1){
        if(strcmp(argv[1], "Adem") == 0){
            algebra = (Algebra*) AdemAlgebra_construct(p, generic, false);
        } else if(strcmp(argv[1], "Milnor") == 0){
            algebra = (Algebra*)MilnorAlgebra_construct(p, generic, NULL);
        } else {
            printf("Unrecognized algebra '%s', using Adem by default.\n", argv[1]);
            algebra = (Algebra*) AdemAlgebra_construct(p, generic, false);
        }
    } else {
        algebra = (Algebra*) AdemAlgebra_construct(p, generic, false);
    }

    uint min_degree = 0;
    uint degree_range = 0;
    uint max_generator_degree = degree_range + min_degree + 1;
    uint graded_dimension[5] = {1};
    Algebra_computeBasis(algebra, degree - min_degree);
    FiniteDimensionalModule *module = 
        FiniteDimensionalModule_construct(algebra, "my_module", min_degree, max_generator_degree, graded_dimension);
    // uint output[1] = {1};
    // FiniteDimensionalModule_setAction(module, 1, 0, min_degree, 1, output);
    // FiniteDimensionalModule_setAction(module, 1, 0, min_degree, 1, output);
    // degree--;
    Resolution *res = Resolution_construct(module, degree, NULL, NULL);
    Resolution_resolveThroughDegree(res, degree);
    // SerializedResolution *sres = Resolution_serialize(res);
    // Resolution *res2 = Resolution_deserialize(module, sres, NULL, NULL);
    // Matrix_print(res->differentials[2]->coimage_to_image_isomorphism[9]);
    // printf("source_dim: %d, columns: %d\n", 
    //     Module_getDimension(&res->differentials[2]->source->module, 9), 
    //     res->differentials[2]->coimage_to_image_isomorphism[9]->columns
    // );

    uint input_homological_degree = 2;
    uint input_degree = 12;
    ResolutionHomomorphism *f = ResolutionHomomorphism_construct(res, res, input_homological_degree, input_degree);
    Vector *v = Vector_construct(p, 1, 0);
    Vector_setEntry(v, 0, 1);
    ResolutionHomomorphism_setBaseMap(f, input_degree, 0, v);
    ResolutionHomomorphism_baseMapReady(f, 1000);

    ResolutionHomomorphism_extend(f, 6, 15);

    // uint idx = FreeModule_operationGeneratorToIndex(M, 0, 0, 12, 0);
    int hom_deg = 2;
    int gen_deg = 4;
    int gen_idx = 0;
    int out_deg = gen_deg - f->internal_degree_shift;
    FreeModule *M = res->modules[hom_deg];
    Vector *output = Vector_construct(p, Module_getDimension((Module*)M,  out_deg), 0);
    // Vector_setEntry(output, idx, 1);
    FreeModuleHomomorphism_applyToGenerator(f->maps[hom_deg - 1], output, 1, gen_deg, gen_idx);
    char buffer[1000];
    FreeModule_element_toJSONString(buffer, (FreeModule*)f->maps[hom_deg - 1]->target,out_deg, output);
    printf("(%d, %d) ==> %s\n", hom_deg, gen_deg, buffer);
    // Matrix_printSlice(res2->differentials[2]->coimage_to_image_isomorphism[8],15,64);
    
    // for(int i = degree - 1 - res->min_degree; i >= 0; i--){
    //     printf("stage %*d: ", 2, i);
    //     array_print("%s\n", &res->modules[i+1]->number_of_generators_in_degree[i], degree - i - res->min_degree);
    // }       
    // MilnorAlgebra_free((MilnorAlgebra*)res->algebra);
    Resolution_free(res);
    FiniteDimensionalModule_free(module);
    return 0;   
}
//*/

