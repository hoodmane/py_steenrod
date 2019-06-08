import os,sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from ctypes_wrap import *
import cFpVector
import cMilnorAlgebra
import steenrod
import steenrod_module
from FreeModule import *



def toC(module):
    module.generate_milnor_action()
    max_degree = max(module.gens.values())
    module.c_algebra = cMilnorAlgebra.cMilnorAlgebra(p=module.milnor_algebra.p, max_degree=max_degree+10)
    number_of_basis_elements_in_degree = [0] * (max_degree + 1)
    basis_element_indices = {}
    index_to_basis_element = {}
    for (b, degree) in module.gens.items():
        index = number_of_basis_elements_in_degree[degree]
        number_of_basis_elements_in_degree[degree] += 1
        basis_element_indices[b] = index
        index_to_basis_element[(degree, index)] = b
    c_list_of_uints_type = (max_degree + 1) * c_uint
    c_number_of_basis_elements_in_degree = c_list_of_uints_type(*number_of_basis_elements_in_degree)

    module.c_algebra.generateBasis(max_degree)
    c_module = CSteenrod.FiniteDimensionalModule_construct(module.c_algebra.c_alg_ptr, max_degree, c_number_of_basis_elements_in_degree)
    module.c_module = c_module
    module.basis_element_indices = basis_element_indices
    module.index_to_basis_element = index_to_basis_element
    module.number_of_basis_elements_in_degree = number_of_basis_elements_in_degree
    for ((op, input), output) in module.milnor_actions.items():
        input_degree = module.gens[input]
        input_index = basis_element_indices[input]
        output_degree = output.degree()
        output_vector = [None] * number_of_basis_elements_in_degree[output_degree]
        op_degree = output_degree - input_degree
        op_index = cMilnorBasisElement.toIndex(algebra, op)        
        for (b, coeff) in output.items():
            output_vector[basis_element_indices[b]] = coeff
        c_vector = cFpVector.cVector(module.p, vector=output_vector)
        CSteenrod.FiniteDimensionalModule_setAction(
            c_module, 
            op_degree, op_index,
            input_degree, input_index,
            c_vector
        )
        c_vector.free()
    return module

def act(module, op, gen):
    c_module = module.c_module
    gen_deg = gen.degree()
    gen_idx = module.basis_element_indices[gen]
    op_deg = op.degree()
    output_degree = op_deg + gen_deg    
    if(output_degree > len(module.number_of_basis_elements_in_degree)):
        return 0
    output_dimension = module.number_of_basis_elements_in_degree[output_degree]
    c_result = cFpVector.cVector(module.p, output_dimension)
    for (op_basis_elt, c) in op.items():
        op_idx = cmilnor.milnor_basis_elt_to_C_index(module.milnor_algebra, op_basis_elt)
        CSteenrod.FiniteDimensionalModule_actOnBasis(cast(c_module, POINTER(c_Module)), c_result.vector, c, op_deg, op_idx, gen_deg, gen_idx)
    result = c_result.unpack()
    c_result.free()
    result = {module.index_to_basis_element[(output_degree, i)] : c for i, c in enumerate(result) if c != 0}
    result = module.get_element(result)
    return result