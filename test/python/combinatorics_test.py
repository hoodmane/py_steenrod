import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,"python"))

from combinatorics import *

def test_base_p_expansion():
    assert base_p_expansion(8,3) == [2,2]
    assert base_p_expansion(33,5) == [3, 1, 1]

def test_binomial():
    assert direct_binomial(21, 2, 23) == 210 % 23
    assert direct_binomial(13, 9, 23) == 715 % 23
    assert direct_binomial(12, 8, 23) == 495 % 23
    assert direct_binomial(13, 8, 23) == 1287 % 23
    assert direct_binomial(14, 8, 23) == 3003 % 23
    assert direct_binomial(14, 9, 23) == 2002 % 23
    assert direct_binomial(15, 5, 23) == 3003 % 23 
    assert direct_binomial(15, 8, 23) == 6435 % 23
    assert direct_binomial(15, 9, 23) == 5005 % 23
    assert direct_binomial(16, 9, 23) == 11440 % 23


def test_multinomial_2():
    assert multinomial_2([1, 2]) == 1
    assert multinomial_2([1, 3]) == 0
    assert multinomial_2([1, 4]) == 1
    assert multinomial_2([2, 4]) == 1
    assert multinomial_2([1, 5]) == 0
    assert multinomial_2([2, 5]) == 1
    assert multinomial_2([2, 6]) == 0
    assert multinomial_2([2, 4, 8]) == 1
    
    
def test_binomial_2():
    assert binomial_2(4, 2) == 0
    assert binomial_2(72, 46) == 0
    assert binomial_2(82, 66) == 1
    assert binomial_2(165, 132) == 1
    assert binomial_2(169, 140) == 0


def test_multinomial_odd():
    assert multinomial_odd([1090, 730], 3) == 1
    assert multinomial_odd([108054, 758], 23) == 18
    assert multinomial_odd([3, 2], 7) == 3

def test_binomial_p():
    pass

def test_xi_degrees():
    assert xi_degrees(17, p=2) == [1, 3, 7, 15]
    assert xi_degrees(17, p=3) == [1, 4, 13]
    assert xi_degrees(400, p=17) == [1, 18, 307]

def test_q_degrees():
    assert tau_degrees(17, p=2) == [1, 3, 7, 15]
    assert tau_degrees(17, p=3) == [1, 5, 17]
    assert tau_degrees(600, p=17) == [1, 33, 577]


def to_set_of_tuples(l):
    return set(map(tuple, l))

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
    assert to_set_of_tuples(WeightedIntegerVectors(10, [1, 4]))   == to_set_of_tuples([[10, 0], [6, 1], [2, 2]])
    assert to_set_of_tuples(WeightedIntegerVectors(7, [1, 3, 7])) == to_set_of_tuples([[7, 0, 0], [4, 1, 0], [1, 2, 0], [0, 0, 1]])
    assert to_set_of_tuples(WeightedIntegerVectors(20, [1, 4, 13])) == to_set_of_tuples([[20, 0, 0], [16, 1, 0], [12, 2, 0], [8, 3, 0], [7, 0, 1], [4, 4, 0], [3, 1, 1], [0, 5, 0]])



def test_restricted_partitions():
    assert to_set_of_tuples(restricted_partitions(8, [1])) == to_set_of_tuples([])
    assert to_set_of_tuples(restricted_partitions(10, [6,4,2])) == to_set_of_tuples([[1,1,0]])
    assert to_set_of_tuples(restricted_partitions(10, [6,4,2,2])) == to_set_of_tuples([[1,1,0,0], [1, 0, 1, 1]])
    
