#include <stdio.h>

#include "main.h"
#include "combinatorics.h"
#include "Algebra.h"
#include "AdemAlgebra.h"
#include "MilnorAlgebra.h"
#include "ResolutionHomomorphism.h"
#include "Resolution.h"


void resolveThroughDegree(Resolution *res, int degree){
    for(int int_deg = res->min_internal_degree; int_deg < degree; int_deg ++){
        for(int hom_deg = 0; hom_deg <= int_deg - res->min_internal_degree; hom_deg++){           
            // printf("(%d, %d)\n", hom_deg, int_deg);
            stepResolution(res, hom_deg, int_deg);
        }
    }
}

void stepResolution(Resolution *res, uint homological_degree, int degree){
    // printf("degree: %d, homological_degree: %d\n", degree, homological_degree);
    // Construct kernel -- say that it's everything.
    // We put the module itself in degree zero and we'll want to hit the whole thing.
    Resolution_step(res, homological_degree, degree);
    // Report the answers.
    // Classes:
    uint num_gens = res->modules[homological_degree + 1]->number_of_generators_in_degree[degree-res->min_internal_degree];
    for(uint i=0; i < num_gens; i++){
        // printf("addClass(%d, %d)\n", homological_degree, degree);
        res->addClass(homological_degree, degree, "");
        if(homological_degree > 0){
            Resolution_computeFiltrationOneProducts(res, homological_degree, degree, i);
        }
    }

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
    resolveThroughDegree(res, degree);
    // SerializedResolution *sres = Resolution_serialize(res);
    // Resolution *res2 = Resolution_deserialize(module, sres, NULL, NULL);
    // Matrix_print(res->differentials[2]->coimage_to_image_isomorphism[9]);
    // printf("source_dim: %d, columns: %d\n", 
    //     Module_getDimension(&res->differentials[2]->source->module, 9), 
    //     res->differentials[2]->coimage_to_image_isomorphism[9]->columns
    // );

    uint input_homological_degree = 1;
    uint input_degree = 2;
    ResolutionHomomorphism *f = ResolutionHomomorphism_construct(res, res, input_homological_degree, input_degree);
    Vector *v = Vector_construct(p, 1, 0);
    Vector_setEntry(v, 0, 1);
    ResolutionHomomorphism_setBaseMap(f, input_degree, 0, v);
    ResolutionHomomorphism_baseMapReady(f, 1000);

    ResolutionHomomorphism_extend(f, 6, 15);

    int hom_deg = 2;
    int gen_deg = 10;
    int gen_idx = 0;
    int out_deg = gen_deg - f->internal_degree_shift;
    FreeModule *M = res->modules[hom_deg];
    Vector *output = Vector_construct(p, Module_getDimension((Module*)M,  out_deg), 0);
    // Vector_setEntry(output, idx, 1);
    FreeModuleHomomorphism_applyToGenerator(f->maps[hom_deg - 1], output, 1, gen_deg, gen_idx);
    char buffer[1000];
    FreeModule_element_toJSONString(buffer, (FreeModule*)f->maps[hom_deg - 1]->target,out_deg, output);
    printf("(%d, %d) ==> %s\n", hom_deg, gen_deg, buffer);
    
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

