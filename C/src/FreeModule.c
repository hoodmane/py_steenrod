#include <assert.h>
#include <stdio.h>

#include "FreeModule.h"

typedef struct {
    Module module;
    uint *number_of_generators_in_degree;
    // private fields
    uint computed_degree;
    FreeModuleOperationGeneratorPair **basis_element_to_opgen_table;
    uint ***generator_to_index_table;
} FreeModuleInternal;

FreeModule *FreeModule_construct(Algebra *algebra, int min_degree, int max_degree){
    uint num_degrees = max_degree - min_degree;
    size_t module_size =  sizeof(FreeModuleInternal) 
        + num_degrees * sizeof(uint) // number_of_generators_in_degree
        + num_degrees * sizeof(FreeModuleOperationGeneratorPair*) //basis_element_to_opgen_table
        + num_degrees * sizeof(uint **) // generator_to_index_table
    ;
    FreeModuleInternal *module = malloc(module_size);
    memset(module, 0, module_size);
    module->module.algebra = algebra;
    module->module.p = algebra->p;
    module->module.computeBasis = FreeModule_computeBasis;
    module->module.getDimension = FreeModule_getDimension;
    module->module.actOnBasis = FreeModule_actOnBasis;
    module->module.max_degree = max_degree;
    module->module.min_degree = min_degree;
    module->computed_degree = 0;
    module->number_of_generators_in_degree = (uint *)(module + 1);
    module->basis_element_to_opgen_table = (FreeModuleOperationGeneratorPair**)(
        module->number_of_generators_in_degree + num_degrees
    );    
    module->generator_to_index_table = (uint ***)(module->basis_element_to_opgen_table + num_degrees);
    assert((char*)(module->generator_to_index_table + num_degrees) == (char*)module + module_size);
    return (FreeModule*)module;
}

void FreeModule_free(FreeModule *module){
    if(module == NULL){
        return;
    }
    FreeModuleInternal *M = (FreeModuleInternal*) module;
    for(uint i = 0; i < module->module.max_degree - module->module.min_degree; i++){
        free(M->generator_to_index_table[i]);
    }
    free(module);
}

bool FreeModule_computeBasis(Module *this __attribute__((unused)), int degree __attribute__((unused))){
    return true;
}

uint FreeModule_getDimension(Module *this, int degree){
    assert(degree >= this->min_degree);
    assert(degree < this->max_degree);
    FreeModule *module = (FreeModule*) this;
    uint result = 0;
    for(int i = this->min_degree; i <= degree; i++){
        // for the excess make sure to use math degree
        result += FreeModule_getNumberOfGensInDegree(module, i)
            * Algebra_getDimension(this->algebra, degree - i, i);
    }
    return result;
}

// Run FreeModule_ConstructBlockOffsetTable(module, degree) before using this on an input in that degree
void FreeModule_actOnBasis(Module *this, Vector *result, uint coeff, int op_deg, uint op_idx, int module_degree, uint module_idx){
    assert(op_idx < Algebra_getDimension(this->algebra, op_deg, module_degree));
    assert(FreeModule_getDimension(this, op_deg + module_degree) <= result->dimension);
    FreeModuleInternal *module = (FreeModuleInternal *) this;
    FreeModuleOperationGeneratorPair operation_generator = FreeModule_indexToOpGen((FreeModule*)module, module_degree, module_idx);
    int module_operation_degree = operation_generator.operation_degree;
    uint module_operation_index = operation_generator.operation_index;
    int generator_degree = operation_generator.generator_degree; 
    uint generator_index  = operation_generator.generator_index;
    // Now all of the output elements are going to be of the form s * x. Find where such things go in the output vector.
    uint num_ops = Algebra_getDimension(this->algebra, module_operation_degree + op_deg, generator_degree);
    uint output_block_min = FreeModule_operationGeneratorToIndex((FreeModule*)module, module_operation_degree + op_deg, 0, generator_degree, generator_index);

    uint output_block_max = output_block_min + num_ops;
    char output_block_memory[Vector_getSize(this->p, 0, 0)];    
    Vector *output_block = Vector_initialize(this->p, output_block_memory, 0, 0);     
    Vector_slice(output_block, result, output_block_min, output_block_max); 
    // Now we multiply s * r and write the result to the appropriate position.
    Algebra_multiplyBasisElements(module->module.algebra, output_block, coeff, op_deg, op_idx, module_operation_degree, module_operation_index, generator_degree);
}


// Compute tables:
//    basis element index     --> operator, generator pair
//    a generator  --> where does that generator's block start?
void FreeModule_ConstructBlockOffsetTable(FreeModule *M, int degree){
    assert(degree>=M->module.min_degree);
    assert(degree < M->module.max_degree);
    int shifted_degree = degree - M->module.min_degree;
    FreeModuleInternal *module = (FreeModuleInternal *) M;
    if(module->generator_to_index_table[shifted_degree] != NULL){
        return;
    }
    // gen_to_idx goes gen_degree => gen_idx => start of block.
    // so gen_to_idx_size should be (number of possible degrees + 1) * sizeof(uint*) + number of gens * sizeof(uint).
    // The other part of the table goes idx => opgen
    // The size should be (number of basis elements in current degree) * sizeof(FreeModuleOperationGeneratorPair)
    // A basis element in degree n comes from a generator in degree i paired with an operation in degree n - i.
    size_t gen_to_idx_size = (degree - M->module.min_degree + 1) * sizeof(uint*);
    size_t total_size = 0;
    for(int gen_deg = M->module.min_degree; gen_deg <= degree; gen_deg++){
        int op_deg = degree - gen_deg;
        uint num_ops = Algebra_getDimension(module->module.algebra, op_deg, gen_deg);
        uint num_gens = FreeModule_getNumberOfGensInDegree(M, gen_deg);
        gen_to_idx_size += num_gens * sizeof(uint);
        total_size += num_gens * num_ops * sizeof(FreeModuleOperationGeneratorPair);
    }
    total_size += gen_to_idx_size;
    char *memory = realloc(module->generator_to_index_table[shifted_degree], total_size);
    uint **generator_to_index_table = (uint**) memory;
    FreeModuleOperationGeneratorPair * basis_element_to_opgen_table = (FreeModuleOperationGeneratorPair *)(memory + gen_to_idx_size);
    module->generator_to_index_table[shifted_degree] = generator_to_index_table;
    module->basis_element_to_opgen_table[shifted_degree] = basis_element_to_opgen_table;
    uint **gentoidx_degree_ptr = generator_to_index_table;
    uint *gentoidx_index_ptr = (uint*)(generator_to_index_table + (shifted_degree + 1));
    uint offset = 0;
    uint generator = 0;
    for(int gen_deg = M->module.min_degree; gen_deg <= degree; gen_deg++){
        *gentoidx_degree_ptr = gentoidx_index_ptr;
        int op_deg = degree - gen_deg;
        uint num_ops = Algebra_getDimension(module->module.algebra, op_deg, gen_deg);
        for(uint gen_idx = 0; gen_idx < FreeModule_getNumberOfGensInDegree(M, gen_deg); gen_idx++){
            *gentoidx_index_ptr = offset;
            for(uint op_idx = 0; op_idx < num_ops; op_idx++){
                basis_element_to_opgen_table->generator_degree = gen_deg;
                basis_element_to_opgen_table->generator_index  = gen_idx;
                basis_element_to_opgen_table->operation_degree = op_deg;
                basis_element_to_opgen_table->operation_index  = op_idx;
                basis_element_to_opgen_table ++;
            }
            gentoidx_index_ptr++;
            generator ++;
            offset += num_ops;
        }
        gentoidx_degree_ptr++;
    }
    assert(gentoidx_degree_ptr == generator_to_index_table + (shifted_degree + 1));
    assert((char*)gentoidx_index_ptr == memory + gen_to_idx_size);
    assert((char*)basis_element_to_opgen_table == memory + total_size);
}



uint FreeModule_getNumberOfGensInDegree(FreeModule *this, int degree){
    return this->number_of_generators_in_degree[degree - this->module.min_degree];
}

uint FreeModule_operationGeneratorToIndex(FreeModule *this, int op_deg, uint op_idx, int gen_deg, uint gen_idx){
    FreeModuleInternal *module = (FreeModuleInternal *)this;
    gen_deg -= this->module.min_degree;
    assert(gen_deg >= 0);
    int deg = op_deg + gen_deg;
    assert(module->generator_to_index_table[deg]!=NULL);
    uint block_idx = module->generator_to_index_table[deg][gen_deg][gen_idx];
    return block_idx + op_idx;
}

FreeModuleOperationGeneratorPair FreeModule_indexToOpGen(FreeModule *this, int degree, uint index){
    degree -= this->module.min_degree;
    FreeModuleInternal *module = (FreeModuleInternal*) this;
    return module->basis_element_to_opgen_table[degree][index];
}

uint FreeModule_element_toJSONString(char *result, FreeModule *this, int degree, Vector *element){
    uint len = 0;
    len += sprintf(result + len, "[");
    for(
        VectorIterator it = Vector_getIterator(element);
        it.has_more;
        it = Vector_stepIterator(it)
    ){
        if(it.value == 0){
            continue;
        }
        len += sprintf(result + len, "{");
        FreeModuleOperationGeneratorPair opgen = FreeModule_indexToOpGen(this, degree, it.index);
        len += sprintf(result + len, "\"op_deg\" : %d,", opgen.operation_degree);
        len += sprintf(result + len, "\"op_idx\" : %d,", opgen.operation_index);
        len += sprintf(result + len, "\"gen_deg\": %d,", opgen.generator_degree);
        len += sprintf(result + len, "\"gen_idx\": %d,", opgen.generator_index);
        len += sprintf(result + len, "\"coeff\" : %d,", it.value);
        len += sprintf(result + len, "\"op_str\" : \"");
        len += Algebra_basisElementToString(this->module.algebra, result + len, opgen.operation_degree, opgen.operation_index);
        len += sprintf(result + len, "\"},");
    }
    len --;
    len += sprintf(result + len, "]");
    return len;
}
