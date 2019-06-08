from ctypes_wrap import *
from cFpVector import *


def toC(algebra, b):
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

def toIndex(algebra, b):
    c_MBE = toC(algebra, b)
    return CSteenrod.MilnorAlgebra_basisElement_toIndex(algebra.c_algebra, c_MBE)

def fromC(algebra, b):
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

def fromIndex(algebra, degree, idx):
    b = CSteenrod.MilnorAlgebra_basisElement_fromIndex(algebra.c_algebra, degree, idx)
    return fromC(algebra, b)