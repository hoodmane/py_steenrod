import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,os.pardir,"python"))
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,os.pardir,"python","CWrappers"))

import pytest
from cAdemAlgebra import *


def check_C_product(a, b):
    return C_product(a,b) == a*b

A2 = cAdemAlgebra(p=2, max_degree=100)
A2gen = cAdemAlgebra(p=2, generic=True, max_degree=20) #100)
A3 = cAdemAlgebra(p=3, max_degree=40)#200)
A5 = cAdemAlgebra(p=5, max_degree=80)#400)

# def test_C_basis():
#     for alg in (A2, A2gen, A3):
#         for dim in range(10):
#             c_basis = [e.toPy() for e in alg.getBasis(dim)]
#             py_basis = alg.py_algebra.basis(dim)
#             assert c_basis == py_basis

product_tests =  [
    (A2, A2.Sq(1), A2.Sq(1)),
    (A2, A2.Sq(2), A2.Sq(2)),
    # (A2, A2.Sq(20,1), A2.Sq(26)),
    (A2, A2.Sq(8,4,1), A2.Sq(6)),
    (A2, A2.Sq(10,4), A2.Sq(4)),
    (A2, A2.Sq(6,1), A2.Sq(8,4))

    # (A2, A2.Sq(1,8), A2.Sq(1,1)),
    # (A2, A2.Sq(6,2), A2.Sq(8,1)),
#
#     (A2gen,A2gen.P(1), A2gen.P(1)),
#     (A2gen,A2gen.P(2), A2gen.P(1)),
#     (A2gen,A2gen.Q(0)*A2gen.P(4), A2gen.Q(0,1)),
# #
#     (A3, A3.P(1),A3.Q(0)),
#     (A3, A3.P(3,3), A3.Q(0,1)),
#     (A3, A3.P(1,1,1), A3.P(1,1,1)),
#     (A3, A3.P(3,3,1), A3.P(1,1,1)),
#     (A3, A3.Q(1) * A3.P(2, 1), A3.Q(0)),
#     (A3, A3.Q(2) * A3.P(1),    A3.Q(0)),
#     (A3, A3.Q(2) * A3.P(3),    A3.Q(0) * A3.P(3)),
#     (A3, A3.Q(1) * A3.P(0, 1), A3.Q(0) * A3.P(3)),
#     (A3, A3.Q(1) * A3.P(1, 1), A3.Q(0) * A3.P(1, 1)),
#     (A3, A3.P(10), A3.P(9)),

#     (A5, A5.Q(1) * A5.P(8, 2),  A5.Q(0) * A5.P(7)),
#     (A5, A5.Q(2) * A5.P(11),  A5.P(4, 2))
]

@pytest.mark.parametrize("A,x,y", product_tests)
def test_C_products(A, x,y):
    c_prod = A.multiply(x, y)
    py_prod = x * y
    assert c_prod == py_prod


def get_random_C_products(num_tests):
    import random    
    algebras = ( A2gen, A3, A5) #A2,
    result = []
    for i in range(num_tests):
        algebra = random.choice(algebras)
        x_deg = 10000; y_deg = 0
        x_dim = 0; y_dim = 0
        max_degree = algebra.max_degree
        while (x_deg + y_deg >= max_degree) or (x_dim == 0) or (y_dim == 0):
            x_deg = random.randint(0,max_degree)
            y_deg = random.randint(0,max_degree)
            x_dim = algebra.getDimension(x_deg)
            y_dim = algebra.getDimension(y_deg)

        x_basis = algebra.py_algebra.basis(x_deg)
        y_basis = algebra.py_algebra.basis(y_deg)
        x = random.choice(x_basis)
        y = random.choice(y_basis)
        result.append((algebra, x, y))
    return result

@pytest.mark.parametrize("A,x,y", get_random_C_products(10))
def test_random_C_products(A,x, y):
    print(A.py_algebra, x, y)
    py_prod = x * y
    prod = A.multiply(x, y)
    assert py_prod == prod

