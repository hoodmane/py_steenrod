import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,"src"))

from combinatorics import *

def test_base_p_expansion():
    pass

def test_binomial():
    pass

def test_multinomial_mod_2():
    pass

def test_binomial_mod2():
    pass


def test_multinomial_odd():
    pass

def test_binomial_modp():
    pass


def test_restricted_partitions():
    assert restricted_partitions(10, [6,4,2]  ) == [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
    assert restricted_partitions(10, [6,4,2,2]) == [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
    assert restricted_partitions(10, [6,4,2] , no_repeats = True) == [[6,4]]
    assert restricted_partitions(10, [6,4,2,2], no_repeats = True) == [[6,4], [6, 2, 2]]

def test_xi_degrees():
    assert xi_degrees(17, p = 2) == [15, 7, 3, 1]
    assert xi_degrees(17, p = 2, reverse = False) == [1, 3, 7, 15]
    assert xi_degrees(17, p = 3) == [13, 4, 1]
    assert xi_degrees(400, p = 17) == [307, 18, 1]


def test_WeightedIntegerVectors():
    """
    Iterate over all ``l`` weighted integer vectors with total weight ``n``.

    INPUT:

    - ``n`` -- an integer
    - ``l`` -- the weights in weakly decreasing order

    EXAMPLES::

        sage: from sage.combinat.integer_vector_weighted import iterator_fast
        sage: list(iterator_fast(3, [2,1,1]))
        [[1, 1, 0], [1, 0, 1], [0, 3, 0], [0, 2, 1], [0, 1, 2], [0, 0, 3]]
        sage: list(iterator_fast(2, [2]))
        [[1]]

    Test that :trac:`20491` is fixed::

        sage: type(list(iterator_fast(2, [2]))[0][0])
        <type 'sage.rings.integer.Integer'&gt
    """
    pass

