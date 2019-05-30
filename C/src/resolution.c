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
    // printf("\n\nhom_deg: %d, mod_deg: %d, num_gens: %d\n", homological_degree, degree, num_gens);
}

void addStructline_doNothing(
    uint source_hom_deg, uint source_int_deg, uint source_idx, 
    uint target_hom_deg, uint target_int_deg, uint target_idx
){
    // printf("\n\nhom_deg: %d, mod_deg: %d, num_gens: %d\n", homological_degree, degree, num_gens);
}


Resolution *Resolution_construct(FiniteDimensionalModule *module, uint max_filtration, uint max_degree){
    Resolution *res = malloc(
        sizeof(Resolution)
        + (max_filtration + 1) * sizeof(FreeModule*)
        + (max_filtration + 1) * sizeof(FreeModuleHomomorphism*)
        + (max_degree + 1) * sizeof(int)
    );
    res->modules = (FreeModule**)(res + 1);
    res->differentials = 
        (FreeModuleHomomorphism**)(res->modules + (max_filtration + 1));
    res->internal_degree_to_resolution_stage = (int*)(res->differentials + max_filtration + 1);
    memset(res->internal_degree_to_resolution_stage, 0, (max_degree+1)*sizeof(int));
    
    res->module = (Module*)module;
    res->algebra= module->module.algebra;
    res->max_degree = max_degree;
    res->modules[0] = (FreeModule*)module;
    res->differentials[0] = FreeModuleHomomorphism_construct((FreeModule*)module, NULL, max_degree);
    res->addClass = addClass_doNothing;
    res->addStructline = addStructline_doNothing;
    return res;
}

void Resolution_resolveThroughDegree(Resolution *res, uint degree){
    for(uint hom_deg = 0; hom_deg < degree; hom_deg++){
       for(uint int_deg = hom_deg; int_deg <= degree; int_deg ++){
            Resolution_step(res, hom_deg, int_deg);
       }
   }
}

void Resolution_step(Resolution *res, uint homological_degree, uint degree){
    // Construct kernel -- say that it's everything.
    // We put the module itself in degree zero and we'll want to hit the whole thing.
    if(homological_degree == 0){
        VectorInterface *vectImpl = &res->algebra->vectorInterface;
        FreeModuleHomomorphism *dminus1 = res->differentials[0];
        uint module_dim = module_getDimension(res->module, degree);
        dminus1->kernel[degree] = Kernel_construct(vectImpl, res->algebra->p, module_dim, module_dim);
        for(uint j = 0; j < module_dim; j++){
            dminus1->kernel[degree]->column_to_pivot_row[j] = j;
            vectImpl->setEntry(dminus1->kernel[degree]->kernel->matrix[j], j, 1);
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
            for(uint j = 0; j < 3; j++){
                uint hj_degree = 1 << j;
                if(hj_degree > degree){
                    break;
                }
                uint gen_degree = degree - hj_degree;
                uint num_target_generators = T->number_of_generators_in_degree[gen_degree];
                for(uint target = 0; target < num_target_generators; target++){
                    uint vector_idx = FreeModule_operationGeneratorToIndex(res->modules[homological_degree], hj_degree, 0, gen_degree, target);
                    if(getVectorEntry(dx, vector_idx) != 0){
                        // There was a product!
                        res->addStructline(homological_degree - 1, gen_degree, target, homological_degree, degree, source);
                    }
                }
            }
        }
    }
}

// Invariants:
void Resolution_generateOldKernelAndComputeNewKernel(Resolution *resolution, uint homological_degree, uint degree){
    VectorInterface *vectImpl = &resolution->algebra->vectorInterface;
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
    uint entries_per_limb = vectImpl->getEntriesPer64Bits(p);
    // Pad the target dimension so that it ends in an aligned position.
    uint padded_target_dimension = 
        ((target_dimension + entries_per_limb - 1)/entries_per_limb)*entries_per_limb;
    uint rows = max(source_dimension, target_dimension);
    uint columns = padded_target_dimension + source_dimension + rows;
    char matrix_memory[getMatrixSize(vectImpl, p, rows, columns)];
    Matrix *full_matrix = initializeMatrix(matrix_memory, vectImpl, p, rows, columns);

    // For the first stage we just want the part of size source_dimension x (padded_target_dimension + source_dimension).
    // Slice matrix out of full_matrix.
    char slice_matrix_memory[Matrix_getSliceSize(source_dimension)];
    Matrix *matrix = Matrix_slice(full_matrix, slice_matrix_memory, 0, source_dimension, 0, padded_target_dimension + source_dimension);
    FreeModuleHomomorphism_getMatrix(current_differential, matrix, degree);
    // Write the identity matrix into the right block
    for(int i = 0; i < source_dimension; i++){
        vectImpl->setEntry(matrix->matrix[i], padded_target_dimension + i, 1);
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
    Kernel *kernel = Kernel_construct(vectImpl, p, kernel_dimension, source_dimension);
    // Write pivots into kernel
    for(uint i = 0; i < source_dimension; i++){
        kernel->column_to_pivot_row[i] = column_to_pivot_row[i + padded_target_dimension] - first_kernel_row;
    }

    // Copy kernel matrix into kernel
    for(uint row = 0; row < kernel_dimension; row++){
        char slice_memory[vectImpl->container_size];
        Vector *slice = (Vector*)slice_memory;
        vectImpl->slice(slice, matrix->matrix[first_kernel_row + row], padded_target_dimension, padded_target_dimension + source_dimension);
        vectImpl->assign(kernel->kernel->matrix[row], slice);
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
            char slice_memory[vectImpl->container_size];
            Vector *slice = (Vector*)slice_memory;
            vectImpl->slice(slice, full_matrix->matrix[current_target_row], 0, previous_cycles->kernel->columns);
            vectImpl->assign(slice, new_image);
            vectImpl->slice(slice, full_matrix->matrix[current_target_row], padded_target_dimension, full_matrix->columns);
            vectImpl->setToZero(slice);
            vectImpl->setEntry(slice, source_dimension + homology_dimension, 1);
            current_target_row++;
            homology_dimension++;
        }
    }
    free(previous_cycles); // This information can now be found in this differential's coimage_to_image_matrix.
    current_differential->source->number_of_generators += homology_dimension;
    current_differential->source->number_of_generators_in_degree[degree] = homology_dimension;
    FreeModuleHomomorphism_AllocateSpaceForNewGenerators(current_differential, degree, homology_dimension);
    for(uint i = 0; i < homology_dimension; i++){
        char slice_memory[vectImpl->container_size]; 
        Vector *slice = (Vector*) slice_memory;
        vectImpl->slice(slice, full_matrix->matrix[first_kernel_row + i], 0, target_dimension);
        FreeModuleHomomorphism_setOutput(current_differential, degree, i, slice);
    }
    FreeModule_ConstructBlockOffsetTable(source, degree);

    // Now the part of the matrix that contains interesting information is current_target_row * (target_dimension + source_dimension + kernel_size).
    // Allocate a matrix coimage_to_image with these dimensions.
    uint coimage_to_image_rows = current_target_row;
    uint coimage_to_image_columns = padded_target_dimension + source_dimension + homology_dimension;
    Matrix *coimage_to_image = constructMatrix(vectImpl, p, coimage_to_image_rows, coimage_to_image_columns);
    current_differential->coimage_to_image_isomorphism[degree] = coimage_to_image;
    // Copy matrix contents to coimage_to_image
    for(uint i = 0; i < coimage_to_image_rows; i++) {
        char slice_memory[vectImpl->container_size]; 
        Vector *slice = (Vector*) slice_memory;
        vectImpl->slice(slice, full_matrix->matrix[i], 0, coimage_to_image_columns);
        vectImpl->assign(coimage_to_image->matrix[i], slice);
    }
    int useless_pivot_row_info[coimage_to_image_columns];
    rowReduce(coimage_to_image, useless_pivot_row_info, 0, 0);
    // printMatrixSlice(coimage_to_image, target_dimension, padded_target_dimension);
    // TODO: assertion about contents of useless_pivot_row_info?
    // Should contain [0,1,2,3,...,n,-1,-1,-1,..., -1].
    resolution->internal_degree_to_resolution_stage[degree] ++;
}

void test2(void (*fn)(uint x)){
    printf("Calling argument...\n");
    fn(43);
    printf("Called it.\n");
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
    MilnorAlgebra *A = constructMilnorAlgebra(2, false, NULL);
    Algebra *algebra = (Algebra*) A;
    algebra_computeBasis(algebra, degree + 10);
    // Vector * v = constructVector2(2, 65, 0);
    // uint i = 64;
    // setVectorEntry(v, i, 1);
    // // setVectorEntry(v, 65, 1);
    // printf("entry%d: %d\n",i, getVectorEntry(v, i));
    // // printf("entry64: %d\n", getVectorEntry(v, 64));
    // printVector(v);
    // printf("\n");
    // return 0;
    uint max_generator_degree = 0;
    uint number_of_generators_in_degree[1] = {1};
    FiniteDimensionalModule *module = 
        FiniteDimensionalModule_construct(algebra, max_generator_degree, number_of_generators_in_degree);
    Resolution *res = Resolution_construct(module, degree, degree);
    if(addClass != NULL){
        res->addClass = addClass;
    }
    if(addStructline != NULL){
        res->addStructline = addStructline;
    }    
    // uint max_deg = 25;
    // uint max_int_deg = 18;  //max_deg;
    // uint max_hom_deg = 3;   //max_deg;
    // for(uint int_deg = 0; int_deg < max_int_deg; int_deg ++){
    //     for(uint hom_deg = 0; hom_deg < max_hom_deg && hom_deg <= int_deg; hom_deg++){
    //         Resolution_step(res, hom_deg, int_deg);
    //     }    
    // }
    uint max_int_deg = degree;
    uint max_hom_deg = degree;
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
    MilnorAlgebra *alg = constructMilnorAlgebra(2, false, NULL);
    GenerateMilnorBasis(alg, 30);
    Vector * result = constructVector2(2, 23, 0);
    MilnorProduct((Algebra*)alg, result, 1, 12, 2, 11, 1);
    // Resolution *res = doResolution(30, NULL, NULL);
    // printf("max_gen_deg: %d\n", res->modules[2]->max_generator_degree);
    // printArray(res->modules[2]->number_of_generators_in_degree, res->modules[2]->max_generator_degree);
    // printf("\n");

    // FreeModuleHomomorphism_applyToBasisElement(res->differentials[7], result, 1, 28, 24);
    printVector(result);
    return 0;
}
//**/

