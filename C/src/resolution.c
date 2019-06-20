//
// Created by Hood on 5/20/2019.
//

// TODO: split up this file a bit.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "combinatorics.h"
#include "algebra.h"
#include "milnor.h"
#include "adem.h"
#include "modules.h"
#include "resolution.h"

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

Resolution * Resolution_construct(
    FiniteDimensionalModule *module, 
    int max_degree,
    void (*addClass)(uint hom_deg, int int_deg, char *cocycle_name),
    void (*addStructline)(
        uint source_hom_deg, int source_int_deg, uint source_idx, 
        uint target_hom_deg, int target_int_deg, uint target_idx
    )    
){
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
    res->differentials[0] = FreeModuleHomomorphism_construct((FreeModule*)module, NULL, max_degree);
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
    for(uint i = 0; i < res->max_degree; i++){
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
        for(uint hom_deg = 0; hom_deg <= int_deg - res->min_degree; hom_deg++){           
            // printf("(%d, %d)\n", hom_deg, int_deg);
            Resolution_step(res, hom_deg, int_deg);
        }
    }
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
        uint module_dim = module_getDimension(res->module, degree);
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
        if((int)op_degree + res->min_degree > degree){
            break;
        }
        printf("hi\n");
        if(op_degree > degree - op_degree){
            printf("hi\n");
            break;
        }

        uint gen_degree = degree - op_degree;
        uint num_target_generators = T->number_of_generators_in_degree[gen_degree - res->min_degree];
        for(uint target = 0; target < num_target_generators; target++){
            uint vector_idx = FreeModule_operationGeneratorToIndex(res->modules[homological_degree], op_degree, op_index, gen_degree, target);
            if(vector_idx >= dx->dimension){
                printf("Out of bounds index when computing product:\n");
                printf("  ==  degree: %d, hom_deg: %d, dim: %d, idx: %d\n", degree, homological_degree, dx->dimension, vector_idx);
            } else {
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
    assert(degree >= (int)homological_degree + resolution->min_degree);
    // printf("degree: %d, resolution->max_degree: %d\n", degree, resolution->max_degree);
    assert(degree < resolution->max_degree);
    uint shifted_degree = degree - resolution->module->min_degree;
    uint p = resolution->algebra->p;
    FreeModuleHomomorphism *current_differential  = resolution->differentials[homological_degree + 1];
    FreeModuleHomomorphism *previous_differential = resolution->differentials[homological_degree];
    FreeModule *source = current_differential->source;

    uint source_dimension = module_getDimension(&current_differential->source->module, degree);
    uint target_dimension = module_getDimension(current_differential->target, degree);
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
    // if(homological_degree <= 2 && degree == 17){
        // printf("\n%d\n",homological_degree);
        // Matrix_printSlice(matrix, target_dimension, padded_target_dimension);
    // }    
    rowReduce(matrix, column_to_pivot_row, 0, 0);//target_dimension, padded_target_dimension);
    // if(homological_degree <= 2 && degree == 17){
        // Matrix_printSlice(matrix, target_dimension, padded_target_dimension);
    // }    


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
        kernel->column_to_pivot_row[i] = column_to_pivot_row[i + padded_target_dimension] - first_kernel_row;
    }

    // Copy kernel matrix into kernel
    for(uint row = 0; row < kernel_dimension; row++){
        char slice_memory[Vector_getContainerSize(p)];
        Vector *slice = Vector_initialize(p, slice_memory, NULL, 0, 0);
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
            char slice_memory[Vector_getContainerSize(p)];
            Vector *slice = Vector_initialize(p, slice_memory, NULL, 0, 0);
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
    Kernel_free(previous_cycles); // This information can now be found in this differential's coimage_to_image_matrix.
    previous_differential->kernel[shifted_degree] = NULL;
    source->number_of_generators_in_degree[shifted_degree] = homology_dimension;
    // Copy the outputs, currently stored in the coimage_to_image matrix, to the FreeModuleHomomorphism outputs field.
    FreeModuleHomomorphism_AllocateSpaceForNewGenerators(current_differential, degree, homology_dimension);
    for(uint i = 0; i < homology_dimension; i++){
        char slice_memory[Vector_getContainerSize(p)]; 
        Vector *slice = Vector_initialize(p, slice_memory, NULL, 0, 0);
        Vector_slice(slice, full_matrix->matrix[first_kernel_row + i], 0, target_dimension);
        FreeModuleHomomorphism_setOutput(current_differential, degree, i, slice);
    }
    FreeModule_ConstructBlockOffsetTable(source, degree);

    // The part of the matrix that contains interesting information is current_target_row * (target_dimension + source_dimension + kernel_size).
    // Allocate a matrix coimage_to_image with these dimensions.
    uint coimage_to_image_rows = current_target_row;
    uint coimage_to_image_columns = padded_target_dimension + source_dimension + homology_dimension;
    Matrix *coimage_to_image = Matrix_construct(p, coimage_to_image_rows, coimage_to_image_columns);
    current_differential->coimage_to_image_isomorphism[shifted_degree] = coimage_to_image;
    // Copy matrix contents to coimage_to_image
    for(uint i = 0; i < coimage_to_image_rows; i++) {
        char slice_memory[Vector_getContainerSize(p)]; 
        Vector *slice = Vector_initialize(p, slice_memory, NULL, 0, 0);
        Vector_slice(slice, full_matrix->matrix[i], 0, coimage_to_image_columns);
        Vector_assign(coimage_to_image->matrix[i], slice);
    }
    // int useless_pivot_row_info[coimage_to_image_columns];
    // rowReduce(coimage_to_image, useless_pivot_row_info, 0, 0);
    
    // printMatrixSlice(coimage_to_image, target_dimension, padded_target_dimension);
    // TODO: assertion about contents of useless_pivot_row_info?
    // Should contain [0,1,2,3,...,n,-1,-1,-1,..., -1].
    resolution->internal_degree_to_resolution_stage[shifted_degree] ++;
}

uint Resolution_gradedDimensionString(char *buffer, Resolution *resolution){
    uint len = 0;
    len += sprintf(buffer + len, "[\n");
    for(uint i = 0; i < resolution->max_degree; i++){
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
            algebra = (Algebra*) AdemAlgebra_construct(p, generic, true);
        } else if(strcmp(argv[1], "Milnor") == 0){
            algebra = (Algebra*)MilnorAlgebra_construct(p, generic, NULL);
        } else {
            printf("Unrecognized algebra '%s', using Adem by default.\n", argv[1]);
            algebra = (Algebra*) AdemAlgebra_construct(p, generic, false);
        }
    } else {
        algebra = (Algebra*) AdemAlgebra_construct(p, generic, false);
    }

    uint min_degree = -3;
    uint degree_range = 1;
    uint max_generator_degree = degree_range + min_degree + 1;
    uint graded_dimension[5] = {2, 1};
    algebra_computeBasis(algebra, degree - min_degree);
    FiniteDimensionalModule *module = 
        FiniteDimensionalModule_construct(algebra, min_degree, max_generator_degree, graded_dimension);
    uint output[1] = {1};
    FiniteDimensionalModule_setAction(module, 1, 0, min_degree, 1, output);
    // FiniteDimensionalModule_setAction(module, 1, 0, min_degree, 1, output);
    // degree--;
    Resolution *res = Resolution_construct(module, degree, NULL, NULL);
    Resolution_resolveThroughDegree(res, degree);
    for(int i = degree - 1 - res->min_degree; i >= 0; i--){
        printf("stage %*d: ", 2, i);
        array_print("%s\n", &res->modules[i+1]->number_of_generators_in_degree[i], degree - i - res->min_degree);
    }       
    res->module = NULL;
    // MilnorAlgebra_free((MilnorAlgebra*)res->algebra);
    res->algebra = NULL;
    Resolution_free(res);
    FiniteDimensionalModule_free(module);
    return 0;
}
//*/

