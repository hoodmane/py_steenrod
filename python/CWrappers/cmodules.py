import os,sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from ctypes_wrap import *
import cFpVector
import cmilnor
import steenrod
import steenrod_module
from FreeModule import *

def FiniteSteenrodModule_to_C(module):
    module.generate_milnor_action()
    max_degree = max(module.gens.values())    
    module.c_algebra = cmilnor.makeCMilnorAlgebra(p=module.milnor_algebra.p, degree=max_degree+10)
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
    module.c_module = c_module
    module.c_algebra = c_algebra
    module.basis_element_indices = basis_element_indices
    module.index_to_basis_element = index_to_basis_element
    module.number_of_basis_elements_in_degree = number_of_basis_elements_in_degree
    for ((op, input), output) in module.milnor_actions.items():
        input_degree = module.gens[input]
        input_index = basis_element_indices[input]
        output_degree = output.degree()
        output_vector = [None] * number_of_basis_elements_in_degree[output_degree]
        for (b, coeff) in output.items():
            output_vector[basis_element_indices[b]] = coeff
        c_vector = cFpVector.cVector(module.p, vector=output_vector)
        op_degree = output_degree - input_degree
        op_index = cmilnor.milnor_basis_elt_to_C_index(algebra, op)
        CSteenrod.FiniteDimensionalModule_setAction(
            c_module, 
            op_degree, op_index,
            input_degree, input_index,
            c_vector
        )
        c_vector.free()
    return module

def c_act_on_fdmodule(module, op, gen):
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
    module.c_module = c_module
    module.c_algebra = cmilnor.makeCMilnorAlgebra(p=module.algebra.p, degree=max_degree+10)    
    module.generator_indices = generator_indices
    module.index_to_generator = index_to_generator
    return module

def FreeModule_from_c(cF, algebra, name_prefix):
    max_generator_degree = cF.contents.max_generator_degree

    number_of_generators_in_degree = [None] * (max_generator_degree + 1)
    for i in range(max_generator_degree + 1):
        number_of_generators_in_degree[i] = cF.contents.number_of_generators_in_degree[i]    
    generator_indices = {}
    index_to_generator = {}

    F = FreeModule(algebra=algebra)
    F.c_module = cF
    F.c_algebra = cmilnor.makeCMilnorAlgebra(p=algebra.p, degree=max_generator_degree+10)        
    F.generator_indices = generator_indices
    F.index_to_generator = index_to_generator

    for i in range(max_generator_degree + 1):
        for j in range(number_of_generators_in_degree[i]):
            gen_name = name_prefix + str(i)
            if number_of_generators_in_degree[i] > 1:
                gen_name += str(j)
            generator_indices[gen_name] = j
            index_to_generator[(i, j)] = gen_name
            F.add_generator(gen_name, i)
    return F

def c_free_module_index_to_py_opgen(module, degree, idx):
    c_module_cast = cast(module.c_module,POINTER(c_Module))
    dimension = CSteenrod.FreeModule_getDimension(c_module_cast, degree)
    if idx >= dimension:
        raise IndexError("Index %s out of bounds. Dimension of module is %s." % (idx, dimension))
    CSteenrod.FreeModule_ConstructBlockOffsetTable(module.c_module, degree)
    c_opgen = CSteenrod.FreeModule_indexToOpGen(module.c_module, degree, idx)
    b = cmilnor.milnor_basis_elt_from_C_idx(module.c_algebra, c_opgen.operation_degree, c_opgen.operation_index)
    gen = module.index_to_generator[(c_opgen.generator_degree, c_opgen.generator_index)]
    return module.get_basis_element(module.algebra.get_basis_element(b), gen)

# def py_basis_elt_to_c_free_module_index(module, elt):
#     c_module_cast = cast(module.c_module,POINTER(c_Module))



def free_module_elt_from_array(module, degree, result_array):
    result = {}
    for (i,c) in enumerate(result_array):
        if c==0:
            continue
        py_opgen = next(iter(c_free_module_index_to_py_opgen(module, degree, i)))
        result[py_opgen] = c
    return module.get_element(result)

def free_module_elt_from_c(module, degree, c_elt):
    result_array = cFpVector.cVector(vector=c_elt).unpack()
    return free_module_elt_from_array(module, degree, result_array)

def c_act_on_free_module(module, op, element):
    c_module = module.c_module
    c_module_cast = cast(c_module, POINTER(c_Module))
    op_deg = op.degree()
    elt_deg = element.degree()
    output_degree = elt_deg + op_deg
    output_dimension = CSteenrod.FreeModule_getDimension(c_module_cast, output_degree)
    c_result = cFpVector.cVector(module.p, output_dimension)
    CSteenrod.FreeModule_ConstructBlockOffsetTable(c_module, elt_deg)
    CSteenrod.FreeModule_ConstructBlockOffsetTable(c_module, output_degree)
    for (op_basis_elt, c1) in op.items():
        op_idx = cmilnor.milnor_basis_elt_to_C_index(module.c_algebra, op_basis_elt)
        for ((elt_op, gen), c2) in element.items():
            elt_op_deg = elt_op.degree()
            elt_op = next(iter(elt_op))
            elt_op_idx = cmilnor.milnor_basis_elt_to_C_index(module.c_algebra, elt_op)
            gen_deg = module.gens[gen]
            gen_idx = module.generator_indices[gen]
            elt_idx = CSteenrod.FreeModule_operationGeneratorToIndex(c_module, elt_op_deg, elt_op_idx, gen_deg, gen_idx)
            CSteenrod.FreeModule_actOnBasis(c_module_cast, c_result.vector, c1*c2, op_deg, op_idx, elt_deg, elt_idx)    
    result = free_module_elt_from_c(module, output_degree, c_result)
    c_result.free()
    return result

def homomorphism_to_c(f, c_S, c_T):
    c_f = CSteenrod.FreeModuleHomomorphism_construct(c_S, cast(c_T, POINTER(c_Module)), 20)
    for i in range(c_S.contents.max_generator_degree+1):
        CSteenrod.FreeModuleHomomorphism_AllocateSpaceForNewGenerators(
            f, i, c_S.contents.number_of_basis_elements_in_degree[i]
        )
    
    for (g, deg) in c_S.py_module.gens.items():
        CSteenrod.FreeModule_getDimension(sour)
        c_output = cFpVector.vector_to_C()
        CSteenrod.FreeModule_setOutput(f, g, c_S.generator_indices[g], )

class cHomomorphism:
    def __init__(self, *, cf, source, target):
        self.cf = cf
        self.source = source
        self.target = target
    
    def __call__(self, input):
        return c_apply_homomorphism(self, input)


def homomorphism_from_c(cf, S, T):
    f = cHomomorphism(cf=cf, source=S, target=T)
    return f

def c_apply_homomorphism(f, element):
    degree = element.degree()
    c_source = f.source.c_module
    c_target_cast = cast(f.target.c_module, POINTER(c_Module))
    output_dimension = CSteenrod.FreeModule_getDimension(c_target_cast, degree)
    c_result = cFpVector.cVector(f.source.p, output_dimension)
    for ((elt_op, gen), coeff) in element.items():
        elt_op_deg = elt_op.degree()
        elt_op = next(iter(elt_op))
        elt_op_idx = cmilnor.milnor_basis_elt_to_C_index(f.source.c_algebra, elt_op)
        gen_deg = f.source.gens[gen]
        gen_idx = f.source.generator_indices[gen]
        elt_idx = CSteenrod.FreeModule_operationGeneratorToIndex(f.source.c_module, elt_op_deg, elt_op_idx, gen_deg, gen_idx)
        # print(elt_op, gen)
        # print("degree: ", degree, "idx: ", elt_idx)
        CSteenrod.FreeModuleHomomorphism_applyToBasisElement(f.cf, c_result.vector, coeff, degree, elt_idx)
        # print(free_module_elt_from_c(f.target, degree, c_result))
    result = free_module_elt_from_c(f.target, degree, c_result.vector)
    c_result.free()
    return result

def c_homomorphism_to_matrix(f, degree):
    c_source_cast = cast(f.source.c_module, POINTER(c_Module))
    c_target_cast = cast(f.target.c_module, POINTER(c_Module))
    input_dimension = CSteenrod.FreeModule_getDimension(c_source_cast, degree)
    output_dimension = CSteenrod.FreeModule_getDimension(c_target_cast, degree)
    cM = cFpVector.cMatrix_construct(2, input_dimension, output_dimension)
    CSteenrod.FreeModuleHomomorphism_getMatrix(f.cf, cM, degree)
    result = cFpVector.matrix_from_C(cM)
    # free matrix
    return result

if __name__ == "__main__":
    A = steenrod.MilnorAlgebra(p=2)
    Sq = A.Sq
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
    c_act_on_fdmodule(cM, Sq(3), x0) # 0
    c_act_on_fdmodule(cM, Sq(0,1), x0) # x3

    T = FreeModule(algebra=A)
    x00 = T.add_generator("x00", 0)
    S = FreeModule(algebra=A)
    x11 = S.add_generator("x11", 1)
    x12 = S.add_generator("x12", 2)
    x14 = S.add_generator("x14", 4)
    x18 = S.add_generator("x18", 8)
    cS = FreeModule_to_c(S, 20)

    d1 = ModuleHomomorphism(S, T)
    d1.add_value("x11", A.Sq(1)*x00)
    d1.add_value("x12", A.Sq(2)*x00)
    d1.add_value("x14", A.Sq(4)*x00)
    d1.add_value("x18", A.Sq(8)*x00)
    # A = F.milnor_algebra