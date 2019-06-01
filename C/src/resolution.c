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
#include "modules.h"
#include "resolution.h"

#define max(a,b) ((a) > (b) ? (a) : (b))
void Resolution_generateOldKernelAndComputeNewKernel(Resolution *resolution, uint homological_degree, uint degree);

void printCallback(Resolution * res, uint homological_degree, uint degree, uint num_gens);

void addClass_doNothing(uint hom_deg, uint int_deg, char *cocycle_name){
    // printf("hom_deg: %d, int_deg: %d\n", hom_deg, int_deg);
}

void addStructline_doNothing(
    uint source_hom_deg, uint source_int_deg, uint source_idx, 
    uint target_hom_deg, uint target_int_deg, uint target_idx
){
    // printf("\n\nhom_deg: %d, mod_deg: %d, num_gens: %d\n", homological_degree, degree, num_gens);
}

Resolution * Resolution_construct(
    FiniteDimensionalModule *module, 
    uint max_degree,
    void (*addClass)(uint hom_deg, uint int_deg, char *cocycle_name),
    void (*addStructline)(
        uint source_hom_deg, uint source_int_deg, uint source_idx, 
        uint target_hom_deg, uint target_int_deg, uint target_idx
    )    
){
    Resolution *res = malloc(
        sizeof(Resolution)
        + (max_degree + 1) * sizeof(FreeModule*)
        + (max_degree + 1) * sizeof(FreeModuleHomomorphism*)
        + (max_degree + 1) * sizeof(int)
    );
    res->modules = (FreeModule**)(res + 1);
    res->differentials = 
        (FreeModuleHomomorphism**)(res->modules + (max_degree + 1));
    res->internal_degree_to_resolution_stage = (int*)(res->differentials + max_degree + 1);
    memset(res->internal_degree_to_resolution_stage, 0, (max_degree+1)*sizeof(int));
    
    res->module = (Module*)module;
    res->algebra= module->module.algebra;
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

void Resolution_resolveThroughDegree(Resolution *res, uint degree){
//     for(uint hom_deg = 0; hom_deg < degree; hom_deg++){
//        for(uint int_deg = hom_deg; int_deg <= degree; int_deg ++){
//             Resolution_step(res, hom_deg, int_deg);
//             printf("(%d, %d)\n", hom_deg, int_deg);
//        }
//    }
    for(uint int_deg = 0; int_deg <= degree; int_deg ++){
        for(uint hom_deg = 0; hom_deg <= int_deg; hom_deg++){           
            // printf("(%d, %d)\n", hom_deg, int_deg);
            Resolution_step(res, hom_deg, int_deg);
       }
   }
}

void Resolution_step(Resolution *res, uint homological_degree, uint degree){
    // Construct kernel -- say that it's everything.
    // We put the module itself in degree zero and we'll want to hit the whole thing.
    if(homological_degree == 0){
        FreeModuleHomomorphism *dminus1 = res->differentials[0];
        uint module_dim = module_getDimension(res->module, degree);
        dminus1->kernel[degree] = Kernel_construct(res->algebra->p, module_dim, module_dim);
        for(uint j = 0; j < module_dim; j++){
            dminus1->kernel[degree]->column_to_pivot_row[j] = j;
            Vector_setEntry(dminus1->kernel[degree]->kernel->matrix[j], j, 1);
        }
    }
    // Construct new FreeModule.
    if(homological_degree == degree){
        res->modules[homological_degree + 1] = 
            FreeModule_construct(res->algebra, res->max_degree, res->max_degree);
        res->differentials[homological_degree + 1] =
            FreeModuleHomomorphism_construct(
                res->modules[homological_degree + 1], 
                (Module*)res->modules[homological_degree],
                res->max_degree
            );
    }
    // Do the work
    Resolution_generateOldKernelAndComputeNewKernel(res, homological_degree, degree);

    // Report the answers.
    // Classes:
    uint num_gens = res->modules[homological_degree + 1]->number_of_generators_in_degree[degree];
    for(uint i=0; i < num_gens; i++){
        res->addClass(homological_degree, degree, "");
    }
    // Products:
    if(homological_degree > 0){
        FreeModuleHomomorphism *d = res->differentials[homological_degree + 1];
        FreeModule *T = (FreeModule*)d->target;        
        for(uint source= 0; source < num_gens; source++){
            Vector *dx = d->outputs[degree][source];
            // 1<<j, j < 3
            for(uint j = 0; j < 3; j++){
                uint hj_degree = 1 << j;
                if(hj_degree > degree){
                    break;
                }
                uint gen_degree = degree - hj_degree;
                uint num_target_generators = T->number_of_generators_in_degree[gen_degree];
                for(uint target = 0; target < num_target_generators; target++){
                    uint vector_idx = FreeModule_operationGeneratorToIndex(res->modules[homological_degree], hj_degree, 0, gen_degree, target);
                    if(vector_idx >= dx->dimension){
                        printf("degree: %d, hom_deg: %d, dim: %d, idx: %d\n", degree, homological_degree, dx->dimension, vector_idx);

                    } else {
                        if(Vector_getEntry(dx, vector_idx) != 0){
                            // There was a product!
                            res->addStructline(homological_degree - 1, gen_degree, target, homological_degree, degree, source);
                        }
                    }
                }
            }
        }
    }
}

// Invariants:
void Resolution_generateOldKernelAndComputeNewKernel(Resolution *resolution, uint homological_degree, uint degree){
    uint p = resolution->algebra->p;
    FreeModuleHomomorphism *current_differential  = resolution->differentials[homological_degree + 1];
    FreeModuleHomomorphism *previous_differential = resolution->differentials[homological_degree];
    FreeModule *source = current_differential->source;
    uint source_dimension = module_getDimension(&current_differential->source->module, degree);
    uint target_dimension = module_getDimension(current_differential->target, degree);

    // assert(source->max_generator_degree == degree);

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
    for(int i = 0; i < source_dimension; i++){
        Vector_setEntry(matrix->matrix[i], padded_target_dimension + i, 1);
    }
    // Row reduce
    int column_to_pivot_row[matrix->columns];
    rowReduce(matrix, column_to_pivot_row, 0, 0);//target_dimension, padded_target_dimension);

    // Stage 1: Find kernel of current differential
    // Locate first kernel row
    uint first_kernel_row = matrix->rows;
    for(uint i = padded_target_dimension; i < matrix->columns; i ++){
        if(column_to_pivot_row[i] != -1){
            first_kernel_row = column_to_pivot_row[i];
            break;
        }
    }
    uint kernel_dimension = matrix->rows - first_kernel_row;
    Kernel *kernel = Kernel_construct(p, kernel_dimension, source_dimension);
    // Write pivots into kernel
    for(uint i = 0; i < source_dimension; i++){
        kernel->column_to_pivot_row[i] = column_to_pivot_row[i + padded_target_dimension] - first_kernel_row;
    }

    // Copy kernel matrix into kernel
    for(uint row = 0; row < kernel_dimension; row++){
        char slice_memory[Vector_getContainerSize(p)];
        Vector *slice = (Vector*)slice_memory;
        Vector_slice(slice, matrix->matrix[first_kernel_row + row], padded_target_dimension, padded_target_dimension + source_dimension);
        Vector_assign(kernel->kernel->matrix[row], slice);
    }
    current_differential->kernel[degree] = kernel;
    // Stage 2: Hit kernel of previous differential. 
    Kernel *previous_cycles = previous_differential->kernel[degree];
    // if(degree == 17){
    //     printf("    previous_cycles:\n");
    //     printMatrix(previous_cycles->kernel);
    // }
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
            char slice_memory[Vector_getContainerSize(p)];
            Vector *slice = (Vector*)slice_memory;
            Vector_slice(slice, full_matrix->matrix[current_target_row], 0, previous_cycles->kernel->columns);
            Vector_assign(slice, new_image);
            Vector_slice(slice, full_matrix->matrix[current_target_row], padded_target_dimension, full_matrix->columns);
            Vector_setToZero(slice);
            Vector_setEntry(slice, source_dimension + homology_dimension, 1);
            current_target_row++;
            homology_dimension++;
        }
    }
    free(previous_cycles); // This information can now be found in this differential's coimage_to_image_matrix.
    current_differential->source->number_of_generators += homology_dimension;
    current_differential->source->number_of_generators_in_degree[degree] = homology_dimension;
    FreeModuleHomomorphism_AllocateSpaceForNewGenerators(current_differential, degree, homology_dimension);
    for(uint i = 0; i < homology_dimension; i++){
        char slice_memory[Vector_getContainerSize(p)]; 
        Vector *slice = (Vector*) slice_memory;
        Vector_slice(slice, full_matrix->matrix[first_kernel_row + i], 0, target_dimension);
        FreeModuleHomomorphism_setOutput(current_differential, degree, i, slice);
    }
    FreeModule_ConstructBlockOffsetTable(source, degree);

    // Now the part of the matrix that contains interesting information is current_target_row * (target_dimension + source_dimension + kernel_size).
    // Allocate a matrix coimage_to_image with these dimensions.
    uint coimage_to_image_rows = current_target_row;
    uint coimage_to_image_columns = padded_target_dimension + source_dimension + homology_dimension;
    Matrix *coimage_to_image = Matrix_construct(p, coimage_to_image_rows, coimage_to_image_columns);
    current_differential->coimage_to_image_isomorphism[degree] = coimage_to_image;
    // Copy matrix contents to coimage_to_image
    for(uint i = 0; i < coimage_to_image_rows; i++) {
        char slice_memory[Vector_getContainerSize(p)]; 
        Vector *slice = (Vector*) slice_memory;
        Vector_slice(slice, full_matrix->matrix[i], 0, coimage_to_image_columns);
        Vector_assign(coimage_to_image->matrix[i], slice);
    }
    int useless_pivot_row_info[coimage_to_image_columns];
    rowReduce(coimage_to_image, useless_pivot_row_info, 0, 0);
    // printMatrixSlice(coimage_to_image, target_dimension, padded_target_dimension);
    // TODO: assertion about contents of useless_pivot_row_info?
    // Should contain [0,1,2,3,...,n,-1,-1,-1,..., -1].
    resolution->internal_degree_to_resolution_stage[degree] ++;
}



#include "milnor.h"
Resolution *doResolution(
    uint degree, 
    void (*addClass)(uint hom_deg, uint int_deg, char *cocycle_name),
    void (*addStructline)(
        uint source_hom_deg, uint source_int_deg, uint source_idx, 
        uint target_hom_deg, uint target_int_deg, uint target_idx
    )
){
    initializePrime(2);
    MilnorAlgebra *A = MilnorAlgebra_construct(2, false, NULL);
    Algebra *algebra = (Algebra*) A;
    algebra_computeBasis(algebra, degree);

    uint max_generator_degree = 2;
    uint graded_dimension[5] = {1,0,1};
    FiniteDimensionalModule *module = 
        FiniteDimensionalModule_construct(algebra, max_generator_degree, graded_dimension);
    uint output[1] = {1};
    FiniteDimensionalModule_setAction(module, 2, 0, 0, 0, output);
    Resolution *res = Resolution_construct(module, degree, addClass, addStructline);

    Resolution_resolveThroughDegree(res, degree);
    // for(int i = max_hom_deg-1; i >= 0; i--){
    //     printf("stage %*d: ", 2, i);
    //     printArray(&res->modules[i+1]->number_of_generators_in_degree[i], max_int_deg - i);
    //     printf("\n");
    // }
    return res;
}

typedef struct {
    uint a;
    uint b;
    uint *c;
} testStruct;

testStruct *test(){
    testStruct *result = malloc(sizeof(testStruct) + 2*sizeof(int));
    result->a = 50;
    result->b = 12;
    result->c = (uint*)(result + 1);
    result->c[0] = 17;
    result->c[1] = 33;
    return result;
}

/**/
int main(){
    // Algebra * A = (Algebra*)MilnorAlgebra_construct(5, true, NULL);
    // MilnorAlgebra_generateBasis(A, 100);
    // printf("constructed\n");
    // algebra_computeBasis(A, 10);
    Resolution *res = doResolution(30, NULL, NULL);
    return 0;
}
//**/

