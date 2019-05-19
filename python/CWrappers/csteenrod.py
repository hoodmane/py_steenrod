import os,sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from ctypes import *
from ctypes_wrap import *

import steenrod

def construct_c_algebra(algebra):
    if hasattr(algebra, "c_algebra"):
        return
    result = c_MilnorAlgebra()
    result.p = algebra.p
    result.generic = algebra.generic
    CSteenrod.constructMilnorAlgebra(result)
    algebra.c_algebra = result

def milnor_basis_elt_to_C(algebra, b):
    bitstring = 0
    p_part = b
    if algebra.generic:
        for x in b[0]:
            bitstring |= 1 << x
        p_part = b[1]
    p_array = (c_ulong * len(p_part))(*p_part)
    result = c_MilnorBasisElement()
    result.q_degree = algebra.basis_q_degree(b)
    result.q_part = bitstring
    result.p_degree = algebra.basis_p_degree(b)
    result.p_length = len(p_part)
    result.p_part = p_array
    return result


def milnor_basis_elt_from_C(b):
    bitstring = b.q_part
    q_part = ()
    i = 0
    while b.q_part != 0:
        if b.q_part & 1 != 0:
            q_part += (i,)
        i += 1
        b.q_part >>= 1
    p_part = [0] * b.p_length
    for i in range(b.p_length):
        p_part[i] = b.p_part[i]
    return (q_part, tuple(p_part))    

def milnor_elt_from_C(algebra, m):
    idx = c_MonomialIndex()
    idx.degree = m.contents.degree
    result = {}
    for i in range(m.contents.algebra_dimension):
        if m.contents.vector[i] == 0:
            continue
        idx.index = i
        b = CSteenrod.GetMilnorBasisElementFromIndex(algebra.c_algebra, idx)
        result[milnor_basis_elt_from_C(b)] = m.contents.vector[i]
    return steenrod.MilnorElement(result, algebra=algebra)

    
def C_product(m1, m2):
    A = m1.algebra
    ret = CSteenrod.allocateMilnorElement(A.c_algebra, m1.degree() + m2.degree())
    b1 = next(iter(m1))
    b2 = next(iter(m2))
    mc = milnor_basis_elt_to_C(A, b1)
    nc = milnor_basis_elt_to_C(A, b2)
    CSteenrod.MilnorProduct(A.c_algebra, ret, mc, nc)
    x = milnor_elt_from_C(A, ret)
    return x
    
def test_C_product():
    A = steenrod.MilnorAlgebra.getInstance(p=3)
    P = A.P
    Q = A.Q
    construct_c_algebra(A)
    CSteenrod.GenerateMilnorBasis(A.c_algebra, 200)
    tests = [
#        (P(1),Q(0)),
#        (P(1,1,1), P(1,1,1)),
#        (P(3,3,1), P(1,1,1)),
        (P(3,3), Q(0,1))
    ]

    for (m1, m2) in tests:
        c_prod = C_product(m1, m2) 
        py_prod = m1 * m2
        if c_prod != py_prod:
            print("Test failed: %s * %s -- " % (m1, m2))
            print("   c_prod : ", c_prod)
            print("   py_prod : ", py_prod)

if __name__ == "__main__":
    A = steenrod.MilnorAlgebra(p=3)
    P = A.P
    Q = A.Q
    construct_c_algebra(A)
    CSteenrod.GenerateMilnorBasis(A.c_algebra, 200)
    m = A.P(2)
    n = A.Q(1)
    mc = milnor_basis_elt_to_C(A, next(iter(m)))
    nc = milnor_basis_elt_to_C(A, next(iter(n)))
    ret = CSteenrod.allocateMilnorElement(A.c_algebra, m.degree() + n.degree())
    CSteenrod.MilnorProduct(A.c_algebra, ret, mc, nc)
    x = milnor_elt_from_C(A, ret)
    
