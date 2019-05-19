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


def milnor_basis_elt_from_C(algebra, b):
    bitstring = b.q_part
    q_part = ()
    i = 0
    bitstring = b.q_part
    while bitstring != 0:
        if bitstring & 1 != 0:
            q_part += (i,)
        i += 1
        bitstring >>= 1
    p_part = [0] * b.p_length
    for i in range(b.p_length):
        p_part[i] = b.p_part[i]
    if algebra.generic:
        b = (q_part, tuple(p_part))
    else:
        b = tuple(p_part)
    return b

def milnor_elt_from_C(algebra, m):
    idx = c_MonomialIndex()
    idx.degree = m.contents.degree
    result = {}
    for i in range(m.contents.algebra_dimension):
        if m.contents.vector[i] == 0:
            continue
        idx.index = i
        b = CSteenrod.GetMilnorBasisElementFromIndex(algebra.c_algebra, idx)
        result[milnor_basis_elt_from_C(algebra, b)] = m.contents.vector[i]
    return algebra.get_element(result)


def C_basis(algebra, dim):
    A = algebra.c_algebra
    c_basisElementList = A.basis_table[dim]
    result = [None] * int(c_basisElementList.length)
    for i in range(c_basisElementList.length):
        b = c_basisElementList.list[i]
        result[i] = algebra.get_basis_element(milnor_basis_elt_from_C(algebra, b))
    return result

def C_product(m1, m2):
    A = m1.algebra
    ret = CSteenrod.allocateMilnorElement(A.c_algebra, m1.degree() + m2.degree())
    b1 = next(iter(m1))
    b2 = next(iter(m2))
    mc = milnor_basis_elt_to_C(A, b1)
    nc = milnor_basis_elt_to_C(A, b2)
    CSteenrod.MilnorProduct(A.c_algebra, ret, mc, nc)
    x = milnor_elt_from_C(A, ret)
    CSteenrod.freeMilnorElement(ret)
    return x

def makeAlgebra(*, p, generic=None, profile = None, dim = 0):
    A = steenrod.MilnorAlgebra.getInstance(p=p, generic=generic, profile=profile)
    construct_c_algebra(A)
    CSteenrod.GenerateMilnorBasis(A.c_algebra, dim)
    return A

def test_C_product():
    A2 = makeAlgebra(p=2, dim=100)
    A2gen = makeAlgebra(p=2, generic=True, dim=100)
    A3 = makeAlgebra(p=3, dim=200)
    tests = [
        (A2.Sq(1), A2.Sq(1)),
        (A2.Sq(2), A2.Sq(2)),
        (A2.Sq(1,8), A2.Sq(1,1)),
        (A2gen.P(1), A2gen.P(1)),
        (A2gen.P(2), A2gen.P(1)),
        (A2gen.Q(0)*A2gen.P(4), A2gen.Q(0,1)),
        (A3.P(1),A3.Q(0)),
        (A3.P(1,1,1), A3.P(1,1,1)),
        (A3.P(3,3,1), A3.P(1,1,1)),
        (A3.P(3,3), A3.Q(0,1)),
    ]

    for (m1, m2) in tests:
        c_prod = C_product(m1, m2) 
        py_prod = m1 * m2
        if c_prod != py_prod:
            print("Test failed: %s * %s -- " % (m1, m2))
            print("   c_prod : ", c_prod)
            print("   py_prod : ", py_prod)

def test_C_basis():
    A2 = makeAlgebra(p=2, dim=50)
    A2gen = makeAlgebra(p=2, generic=True, dim=50)
    A3 = makeAlgebra(p=3, dim=50)
    
    for alg in (A2, A2gen, A3):
        for dim in range(50):
            c_basis = C_basis(alg, dim)
            py_basis = alg.basis(dim)
            py_basis.reverse()
            if c_basis != py_basis:
                print("Discrepency for algebra %s in dimension %s." % (alg, dim))
                print("  c_basis:", c_basis)
                print("  py_basis:", py_basis)
            

if __name__ == "__main__":
    A2 = makeAlgebra(p=2, dim=100)
    A3 = makeAlgebra(p=3, dim=100)
    
