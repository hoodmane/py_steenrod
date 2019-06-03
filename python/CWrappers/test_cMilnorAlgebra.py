import pytest
import cMilnorAlgebra 


def check_C_product(a, b):
    return C_product(a,b) == a*b

A2 = cMilnorAlgebra.construct(p=2, degree=100)
A2gen = cMilnorAlgebra.construct(p=2, generic=True, degree=100)
A3 = cMilnorAlgebra.construct(p=3, degree=200)
A5 = cMilnorAlgebra.construct(p=5, degree=400)

def test_C_basis():
    for alg in (A2, A2gen, A3):
        for dim in range(50):
            c_basis = cMilnorAlgebra.getBasis(alg, dim)
            py_basis = alg.basis(dim)
            py_basis.reverse()
            assert c_basis == py_basis

product_tests =  [
    (A2.Sq(1), A2.Sq(1)),
    (A2.Sq(2), A2.Sq(2)),
    (A2.Sq(1,8), A2.Sq(1,1)),
    (A2.Sq(6,2), A2.Sq(8,1)),
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
    (A3.P(10), A3.P(9)),

    (A5.Q(1) * A5.P(8, 2),  A5.Q(0) * A5.P(7)),
    (A5.Q(2) * A5.P(11),  A5.P(4, 2))
]

@pytest.mark.parametrize("x,y", product_tests)
def test_C_products(x,y):
    c_prod = cMilnorAlgebra.multiply(x, y)
    py_prod = x * y
    assert c_prod == py_prod


def get_random_C_products(num_tests):
    import random    
    algebras = (A2, A2gen, A3, A5)
    result = []
    for i in range(num_tests):
        algebra = random.choice(algebras)
        x_deg = 10000; y_deg = 0
        x_dim = 0; y_dim = 0
        max_degree = algebra.c_max_degree
        while (x_deg + y_deg > max_degree) or (x_dim == 0) or (y_dim == 0):
            x_deg = random.randint(0,max_degree)
            y_deg = random.randint(0,max_degree)
            x_dim = cMilnorAlgebra.getDimension(algebra, x_deg)
            y_dim = cMilnorAlgebra.getDimension(algebra, y_deg)

        x_basis = cMilnorAlgebra.getBasis(algebra, x_deg)
        y_basis = cMilnorAlgebra.getBasis(algebra, y_deg)
        x = random.choice(x_basis)
        y = random.choice(y_basis)
        result.append((x, y))
    return result

# @pytest.mark.parametrize("x,y", get_random_C_products(100))
# def test_random_C_products(x, y):
#     #print("%s ( %s, %s, %s ) : %s * %s" % (algebra, x_deg, y_deg, x_deg + y_deg,  x, y))
#     py_prod = x * y
#     prod = cMilnorAlgebra.multiply(x, y)
#     assert py_prod == prod

