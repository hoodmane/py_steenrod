import os,sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from ctypes_wrap import *
import cFpVector
import cAlgebra
import cFreeModule
import steenrod
import steenrod_module
from FreeModule import *

def construct(p, max_degree):
    milnor_algebra = cAlgebra.cMilnorAlgebra(p=p, max_degree=max_degree+10)
    c_algebra = cast(milnor_algebra.c_algebra, POINTER(c_Algebra))
    M = CSteenrod.FreeModule_construct(c_algebra, max_degree)
    return M

def toC(module, max_degree=None):
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
    c_module = construct(module.p, max_degree)
    c_module.contents.number_of_generators = sum(number_of_generators_in_degree)
    for (i, n) in enumerate(number_of_generators_in_degree):
        c_module.contents.number_of_generators_in_degree[i] = n
    module.c_module = c_module
    module.c_algebra = cAlgebra.cMilnorAlgebra(p=module.algebra.p, max_degree=max_degree+10)    
    module.generator_indices = generator_indices
    module.index_to_generator = index_to_generator
    return module

def fromC(cF, algebra, name_prefix):
    max_generator_degree = cF.contents.module.max_degree
    number_of_generators_in_degree = [None] * (max_generator_degree)
    for i in range(max_generator_degree ):
        number_of_generators_in_degree[i] = cF.contents.number_of_generators_in_degree[i]    
    generator_indices = {}
    index_to_generator = {}
    F = FreeModule(algebra=algebra.py_algebra)
    F.c_module = cF
    F.c_algebra = algebra
    F.generator_indices = generator_indices
    F.index_to_generator = index_to_generator
    for i in range(max_generator_degree):
        for j in range(number_of_generators_in_degree[i]):
            gen_name = name_prefix + str(i)
            if number_of_generators_in_degree[i] > 1:
                gen_name += str(j)
            generator_indices[gen_name] = j
            index_to_generator[(i, j)] = gen_name
            F.add_generator(gen_name, i)
    return F

def indexToPyOpgen(module, degree, idx):
    c_module_cast = cast(module.c_module,POINTER(c_Module))
    dimension = CSteenrod.FreeModule_getDimension(c_module_cast, degree)
    if idx >= dimension:
        raise IndexError("Index %s out of bounds. Dimension of module is %s." % (idx, dimension))
    CSteenrod.FreeModule_ConstructBlockOffsetTable(module.c_module, degree)
    c_opgen = CSteenrod.FreeModule_indexToOpGen(module.c_module, degree, idx)
    b = module.c_algebra.basisEltFromIdx(c_opgen.operation_degree, c_opgen.operation_index)
    b = b.toPy()
    gen = module.index_to_generator[(c_opgen.generator_degree, c_opgen.generator_index)]
    return module.get_basis_element(b, gen)

# def py_basis_elt_to_c_free_module_index(module, elt):
#     c_module_cast = cast(module.c_module,POINTER(c_Module))



def elementFromArray(module, degree, result_array):
    result = {}
    for (i,c) in enumerate(result_array):
        if c==0:
            continue
        py_opgen = next(iter(cFreeModule.indexToPyOpgen(module, degree, i)))
        result[py_opgen] = c
    return module.get_element(result)

def elementFromC(module, degree, c_elt):
    result_array = cFpVector.cVector(module.p, vector=c_elt).unpack()
    return elementFromArray(module, degree, result_array)

def act(module, op, element):
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
        op_idx = module.algebra.basis_type.toIndex(module.c_algebra, op_basis_elt)
        for ((elt_op, gen), c2) in element.items():
            elt_op_deg = elt_op.degree()
            elt_op = next(iter(elt_op))
            elt_op_idx = module.algebra.basis_type.toIndex(module.c_algebra, elt_op)
            gen_deg = module.gens[gen]
            gen_idx = module.generator_indices[gen]
            elt_idx = CSteenrod.FreeModule_operationGeneratorToIndex(c_module, elt_op_deg, elt_op_idx, gen_deg, gen_idx)
            CSteenrod.FreeModule_actOnBasis(c_module_cast, c_result.vector, c1*c2, op_deg, op_idx, elt_deg, elt_idx)    
    result = free_module_elt_from_c(module, output_degree, c_result)
    c_result.free()
    return result