import os,sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from ctypes import *
from ctypes_wrap import *
from cFpVector import *

import steenrod

def construct_C_algebra(algebra):
    if hasattr(algebra, "c_algebra"):
        return
    algebra.c_algebra = CSteenrod.constructMilnorAlgebra(algebra.p, algebra.generic, POINTER(c_Profile)())
    algebra.c_max_degree = 0

def makeCMilnorAlgebra(*, p, generic=None, profile = None, degree = 0):
    algebra = steenrod.MilnorAlgebra.getInstance(p=p, generic=generic, profile=profile)
    construct_C_algebra(algebra)
    c_GenerateMilnorBasis(algebra, degree)
    return algebra

def c_GenerateMilnorBasis(algebra, degree):
    if degree > algebra.c_max_degree:
        CSteenrod.GenerateMilnorBasis(cast(algebra.c_algebra, POINTER(c_Algebra)), degree)
        algebra.c_max_degree = degree

def milnor_basis_elt_to_C(algebra, b):
    bitstring = 0
    p_part = b
    if algebra.generic:
        for x in b[0]:
            bitstring |= 1 << x
        p_part = b[1]
    p_array = (c_uint * len(p_part))(*p_part)
    result = c_MilnorBasisElement()
    result.q_degree = algebra.basis_q_degree(b)
    result.q_part = bitstring
    result.p_degree = algebra.basis_p_degree(b)
    result.p_length = len(p_part)
    result.p_part = p_array
    return result

def milnor_basis_elt_to_C_index(algebra, b):
    c_MBE = milnor_basis_elt_to_C(algebra, b)
    return CSteenrod.GetIndexFromMilnorBasisElement(algebra.c_algebra, c_MBE)

def milnor_basis_elt_from_C_MBE(algebra, b):
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

def milnor_basis_elt_from_C_idx(algebra, degree, idx):
    b = CSteenrod.GetMilnorBasisElementFromIndex(algebra.c_algebra, degree, idx)
    return milnor_basis_elt_from_C_MBE(algebra, b)

def milnor_elt_from_C(algebra, m):
    result = {}
    for i in range(m.dimension):
        entry = c_getVectorEntry(m, i)
        if entry == 0:
            continue
        result[milnor_basis_elt_from_C_idx(algebra, m.degree, b)] = entry
    return algebra.get_element(result)

def C_basis(algebra, degree):
    if degree > algebra.c_max_degree:
        raise Exception("C basis only known through degree %s < %s." % (algebra.c_max_degree, degree))
    A = algebra.c_algebra
    c_basisElementList = CSteenrod.GetMilnorAlgebraBasis(cast(algebra.c_algebra, POINTER(c_Algebra)), degree)
    result = [None] * int(c_basisElementList.length)
    for i in range(c_basisElementList.length):
        b = c_basisElementList.list[i]
        result[i] = algebra.get_basis_element(milnor_basis_elt_from_C(algebra, b))
    return result

def C_dimension(algebra, degree):
    return CSteenrod.GetMilnorAlgebraDimension(cast(algebra.c_algebra, POINTER(c_Algebra)), degree)

def C_product(m1, m2):
    algebra = m1.algebra
    m1_deg = m1.degree()
    m2_deg = m2.degree()
    out_degree = m1_deg + m2_deg
    if out_degree > algebra.c_max_degree:
        raise Exception("C basis only known through degree %s < %s." % (algebra.c_max_degree, out_degree))
    out_dimension = C_dimension(algebra, out_degree)
    ret = construct_c_vector(algebra.p, out_dimension)
    ret.degree = out_degree
    b1 = next(iter(m1))
    b2 = next(iter(m2))
    m1_idx = milnor_basis_elt_to_C_index(algebra, b1)
    m2_idx = milnor_basis_elt_to_C_index(algebra, b2)
    CSteenrod.MilnorProduct(cast(algebra.c_algebra, POINTER(c_Algebra)), ret, 1, m1_deg, m1_idx, m2_deg, m2_idx)
    x = milnor_elt_from_C(algebra, ret)
    CSteenrod.freeVector(ret)
    return x


            

def check_C_product(a, b):
    return C_product(a,b) == a*b

def test_C_product():
    A2 = makeCMilnorAlgebra(p=2, degree=100)
    A2gen = makeCMilnorAlgebra(p=2, generic=True, degree=100)
    A3 = makeCMilnorAlgebra(p=3, degree=200)
    A5 = makeCMilnorAlgebra(p=5, degree=400)
    tests = [
        (A2.Sq(1), A2.Sq(1)),
        (A2.Sq(2), A2.Sq(2)),
        (A2.Sq(1,8), A2.Sq(1,1)),
#
        (A2gen.P(1), A2gen.P(1)),
        (A2gen.P(2), A2gen.P(1)),
        (A2gen.Q(0)*A2gen.P(4), A2gen.Q(0,1)),
#
        (A3.P(1),A3.Q(0)),
        (A3.P(3,3), A3.Q(0,1)),
        (A3.P(1,1,1), A3.P(1,1,1)),
        (A3.P(3,3,1), A3.P(1,1,1)),
        (A3.Q(1) * A3.P(2, 1), A3.Q(0)),
        (A3.Q(2) * A3.P(1),    A3.Q(0)),
        (A3.Q(2) * A3.P(3),    A3.Q(0) * A3.P(3)),
        (A3.Q(1) * A3.P(0, 1), A3.Q(0) * A3.P(3)),
        (A3.Q(1) * A3.P(1, 1), A3.Q(0) * A3.P(1, 1)),

        (A5.Q(1) * A5.P(8, 2),  A5.Q(0) * A5.P(7)),
        (A5.Q(2) * A5.P(11),  A5.P(4, 2))
    ]

    # p = 5
    # Q(1) P(2) * P(11, 0, 1)
    # Q(2) P(11) * P(4, 2)
    # Q(0) P(3, 0, 1) * Q(0) Q(1) P(1)
    # Q(0) Q(1) Q(2) P(6, 1) * P(4)
    #
    # (A5.Q(1) * A5.P(2), A5.P(1, 1))
    # ( 16, 56, 72 ) : P(2) * P(1, 1)
    #  ( 17, 8, 25 ) : Q(0) P(2) * P(1)
    # ( 56, 16, 72 ) : P(1, 1) * P(2)
    # A2.Sq(28,12), A2.Sq(3,2)
    # ( 60, 10, 70 ) : Sq(22, 8, 2) * Sq(1, 3)
    # ( 34, 13, 47 ) : Sq(4, 10) * Sq(10, 1)
    # ( 32, 12, 44 ) : Sq(8, 8) * Sq(9, 1)
    # 42: Sq(12, 8) * Sq(3, 1)
    #  Sq(16, 4) * Sq(8, 1)
    for (m1, m2) in tests:
        print("Test : %s * %s " % (m1, m2))    
        c_prod = C_product(m1, m2)
        py_prod = m1 * m2
        if c_prod != py_prod:
            print("Test failed: p=%s, %s * %s -- " % (m1.algebra.p, m1, m2))
            print("   c_prod : ", c_prod)
            print("   py_prod : ", py_prod)


    algebras = (A2, A2gen, A3, A5)
    import random
    for i in range(5000):
        algebra = random.choice(algebras)
        x_deg = 10000; y_deg = 0
        x_dim = 0; y_dim = 0
        max_degree = 42 # algebra.c_max_degree
        while (x_deg + y_deg > max_degree) or (x_dim == 0) or (y_dim == 0):
            x_deg = random.randint(0,max_degree)
            y_deg = random.randint(0,max_degree)
            x_dim = C_dimension(algebra, x_deg)
            y_dim = C_dimension(algebra, y_deg)

        x_basis = C_basis(algebra, x_deg)
        y_basis = C_basis(algebra, y_deg)
        x = random.choice(x_basis)
        y = random.choice(y_basis)
        print("%s ( %s, %s, %s ) : %s * %s" % (algebra, x_deg, y_deg, x_deg + y_deg,  x, y))
        py_prod = x * y
        prod = C_product(x, y)
        if py_prod != prod:
            print("%s ( %s, %s, %s ) : %s * %s" % (algebra, x_deg, y_deg, x_deg + y_deg,  x, y))
            print("  Test %s failed" % i)



def test_C_basis():
    A2 = makeCMilnorAlgebra(p=2, degree=50)
    A2gen = makeCMilnorAlgebra(p=2, generic=True, degree=50)
    A3 = makeCMilnorAlgebra(p=3, degree=50)

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
    pass
    A2 = makeCMilnorAlgebra(p=2, degree=100)
    #A3 = makeCMilnorAlgebra(p=3, dim=200)
    
