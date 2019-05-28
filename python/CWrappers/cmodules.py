import os,sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from ctypes_wrap import *
import cFpVector
import cmilnor
import steenrod
import steenrod_module

def FiniteSteenrodModule_to_C(module):
    module.generate_milnor_action()
    max_degree = max(module.gens.values())    
    module.c_milnor_algebra = cmilnor.makeCMilnorAlgebra(p=module.milnor_algebra.p, degree=max_degree+10)
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

    algebra = module.milnor_algebra
    cmilnor.construct_C_algebra(algebra)
    c_algebra = cast(module.milnor_algebra.c_algebra, POINTER(c_Algebra))
    cmilnor.c_GenerateMilnorBasis(module.milnor_algebra, max_degree)
    c_module = CSteenrod.FiniteDimensionalModule_construct(c_algebra, max_degree, c_number_of_basis_elements_in_degree)
    c_module.py_module = module
    c_module.c_milnor_algebra = module.c_milnor_algebra
    c_module.basis_element_indices = basis_element_indices
    c_module.index_to_basis_element = index_to_basis_element
    c_module.number_of_basis_elements_in_degree = number_of_basis_elements_in_degree
    for ((op, input), output) in M.milnor_actions.items():
        input_degree = M.gens[input]
        input_index = basis_element_indices[input]
        output_degree = output.degree()
        output_vector = [None] * number_of_basis_elements_in_degree[output_degree]
        for (b, coeff) in output.items():
            output_vector[basis_element_indices[b]] = coeff
        c_vector = cFpVector.vector_to_C(module.p, output_vector)
        op_degree = output_degree - input_degree
        op_index = cmilnor.milnor_basis_elt_to_C_index(algebra, op)
        CSteenrod.FiniteDimensionalModule_setAction(
            c_module, 
            op_degree, op_index,
            input_degree, input_index,
            c_vector
        )
        CSteenrod.freeVector(c_vector)
    return c_module

def c_act_on_fdmodule(c_module, op, gen):
    gen_deg = gen.degree()
    gen_idx = c_module.basis_element_indices[gen]
    op_deg = op.degree()
    op_idx = cmilnor.milnor_basis_elt_to_C_index(c_module.c_milnor_algebra, next(iter(op)))
    output_degree = op_deg + gen_deg
    if(output_degree > len(c_module.number_of_basis_elements_in_degree)):
        return 0
    output_dimension = c_module.number_of_basis_elements_in_degree[output_degree]
    c_result = cFpVector.construct_c_vector(c_module.py_module.p, output_dimension)
    CSteenrod.FiniteDimensionalModule_actOnBasis(cast(c_module, POINTER(c_Module)), c_result, 1, op_deg, op_idx, gen_deg, gen_idx)
    result = cFpVector.vector_from_C(c_result)
    CSteenrod.freeVector(c_result)
    result = {c_module.index_to_basis_element[(output_degree, i)] : c for i, c in enumerate(result) if c != 0}
    result = c_module.py_module.get_element(result)
    return result

def c_FreeModule_construct(p, max_generator_degree, max_degree=None):
    if max_degree == None:
        max_degree = max_generator_degree
    milnor_algebra = cmilnor.makeCMilnorAlgebra(p=p, degree=max_degree+10)
    c_algebra = cast(milnor_algebra.c_algebra, POINTER(c_Algebra))
    M = CSteenrod.FreeModule_construct(c_algebra, max_generator_degree, max_degree)
    return M

def FreeModule_to_c(module, max_degree=None):
    max_generator_degree = max(module.gens.values())  
    if(max_degree == None):
        max_degree = max_generator_degree
    number_of_generators_in_degree = [0] * (max_generator_degree + 1)
    generator_indices = {}
    index_to_generator = {}
    for (b, degree) in module.gens.items():
        index = number_of_generators_in_degree[degree]
        number_of_generators_in_degree[degree] += 1
        generator_indices[b] = index
        index_to_generator[(degree, index)] = b
    c_module = c_FreeModule_construct(module.p, max_generator_degree, max_degree)
    c_module.contents.number_of_generators = sum(number_of_generators_in_degree)
    for (i, n) in enumerate(number_of_generators_in_degree):
        c_module.contents.number_of_generators_in_degree[i] = n
    c_module.py_module = module
    c_module.generator_indices = generator_indices
    c_module.index_to_generator = index_to_generator
    c_module.c_milnor_algebra = cmilnor.makeCMilnorAlgebra(p=module.milnor_algebra.p, degree=100)    
    return c_module

def c_act_on_free_module(module, op, elt_op, element):
    op_deg = op.degree()
    elt_deg = element.degree() + elt_op.degree()
    out_deg = elt_deg + op_deg
    CSteenrod.FreeModule_ConstructBlockOffsetTable(module, elt_deg)
    CSteenrod.FreeModule_ConstructBlockOffsetTable(module, out_deg)
    op_idx = cmilnor.milnor_basis_elt_to_C_index(module.c_milnor_algebra, next(iter(op)))
    (_, gen) = next(iter(element))
    elt_op_deg = elt_op.degree()
    elt_op_idx = cmilnor.milnor_basis_elt_to_C_index(module.c_milnor_algebra, next(iter(elt_op)))
    gen_deg = module.py_module.gens[gen]
    gen_idx = module.generator_indices[gen]
    elt_idx = CSteenrod.FreeModule_operationGeneratorToIndex(module, elt_op_deg, elt_op_idx, gen_deg, gen_idx)
    module_cast = cast(module, POINTER(c_Module))
    output_dimension = CSteenrod.FreeModule_getDimension(module_cast, out_deg)
    c_result = cFpVector.construct_c_vector(module.py_module.p, output_dimension)
    CSteenrod.FreeModule_actOnBasis(module_cast, c_result, 1, op_deg, op_idx, elt_deg, elt_idx)    
    result = cFpVector.vector_from_C(c_result)
    CSteenrod.freeVector(c_result)
    return result


if __name__ == "__main__":
    # M = steenrod_module.FiniteSteenrodModule(p=2)
    # x0 = M.add_basis_element("x0", 0)
    # x1 = M.add_basis_element("x1", 1)
    # x2 = M.add_basis_element("x2", 2)
    # x3 = M.add_basis_element("x3", 3)
    # x4 = M.add_basis_element("x4", 4)
    # M.add_Sq_action(2, x0, x2)
    # M.add_Sq_action(2, "x2", x4)
    # M.add_Sq_action(1, "x0", x1)
    # M.add_Sq_action(2, "x1", x3)
    # M.add_Sq_action(1, "x3", x4)
    # M.add_Sq_action(3, "x1", x4)
    # M.validate()
    # cM = FiniteSteenrodModule_to_C(M)
    # c_act_on_fdmodule(cM, M.milnor_algebra.Sq(3), x0) # 0
    # c_act_on_fdmodule(cM, M.milnor_algebra.Sq(0,1), x0) # x3
    F = steenrod_module.FreeSteenrodModule(p=2)
    x0 = F.add_generator("x0", 0)
    x1 = F.add_generator("x1", 1)
    cF = FreeModule_to_c(F, 20)
    A = F.milnor_algebra