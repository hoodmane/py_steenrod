import os,sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from ctypes_wrap import *
import cmilnor
import steenrod
import steenrod_module

def FiniteSteenrodModule_to_C(module):
    module.generate_milnor_action()
    
    max_degree = max(module.gens.values())
    number_of_basis_elements_in_degree = [0] * (max_degree + 1)
    basis_element_indices = {}
    for (b, degree) in module.gens.items():
        basis_element_indices[b] = number_of_basis_elements_in_degree[degree]
        number_of_basis_elements_in_degree[degree] += 1
    print("py: ", number_of_basis_elements_in_degree)
    c_list_of_uints_type = (max_degree + 1) * c_uint
    c_number_of_basis_elements_in_degree = c_list_of_uints_type(*number_of_basis_elements_in_degree)
    print("num_deg:", [c_number_of_basis_elements_in_degree[i] for i in range(max_degree+1)])

    print("Constructing Algebra")
    algebra = module.milnor_algebra
    cmilnor.construct_C_algebra(algebra)
    c_algebra = cast(module.milnor_algebra.c_algebra, POINTER(c_Algebra))
    print("Generating Milnor Basis")
    cmilnor.c_GenerateMilnorBasis(module.milnor_algebra, max_degree)
    print("Constructing Module")
    c_module = CSteenrod.constructFiniteDimensionalModule(c_algebra, max_degree, c_number_of_basis_elements_in_degree)
    print("Adding actions to Module")
    for ((op, input), output) in M.milnor_actions.items():
        input_degree = M.gens[input]
        input_index = basis_element_indices[input]
        output_degree = output.degree()
        output_vector_type = number_of_basis_elements_in_degree[output_degree] * c_uint
        output_vector = output_vector_type()
        for (b, coeff) in output.items():
            output_vector[basis_element_indices[b]] = coeff
        
        op_degree = output_degree - input_degree
        op_index = cmilnor.milnor_basis_elt_to_C_index(algebra, op)
        print("Add action")
        CSteenrod.addActionToFiniteDimensionalModule(
            c_module, 
            op_degree, op_index,
            input_degree, input_index,
            output_vector
        )
    return c_module
    
    
if __name__ == "__main__":
    M = steenrod_module.FiniteSteenrodModule(p=2)
    x0 = M.add_basis_element("x0", 0)
    x1 = M.add_basis_element("x1", 1)
    x2 = M.add_basis_element("x2", 2)
    x3 = M.add_basis_element("x3", 3)
    x4 = M.add_basis_element("x4", 4)
    M.add_Sq_action(2, x0, x2)
    M.add_Sq_action(2, "x2", x4)
    M.add_Sq_action(1, "x0", x1)
    M.add_Sq_action(2, "x1", x3)
    M.add_Sq_action(1, "x3", x4)
    M.add_Sq_action(3, "x1", x4)
    M.validate()
    cM = FiniteSteenrodModule_to_C(M)
    
