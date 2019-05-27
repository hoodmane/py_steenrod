
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
        c_prod = C_product(m1, m2)
        py_prod = m1 * m2
        #print("Test : %s * %s " % (m1, m2))
        if c_prod != py_prod:
            print("Test failed: %s * %s -- " % (m1, m2))
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
        #print("%s ( %s, %s, %s ) : %s * %s" % (algebra, x_deg, y_deg, x_deg + y_deg,  x, y))
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
