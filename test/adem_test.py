import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,"src"))

import math
from functools import reduce

from adem import *

#class MinimalAdemAlgebra:
    # ...
    
    
def test_adem_basis_elt_generic_map():
    """
        *, P_fn, b, basis_elt
    """
    P_fn =  lambda P : "P%s" % P
    result = adem_basis_elt_generic_map( P_fn = P_fn, b = "b", basis_elt = [1] )
    assert result == ["b"]
    result = adem_basis_elt_generic_map( P_fn = P_fn, b = "b", basis_elt = [1, 2, 0] )
    assert result == ["b", "P2"]
    result = adem_basis_elt_generic_map( P_fn = P_fn, b = "b", basis_elt = [1, 4, 1, 1, 0] )
    assert result == ["b", "P4", "b", "P1"]
    
def test_adem_basis_elt_2_map_test():
    """
        This function is a little pointless but...
    """
    result = adem_basis_elt_2_map(Sq_fn = lambda Sq : "Sq%s" % Sq, basis_elt = [4, 2, 1])
    assert result == ["Sq4", "Sq2", "Sq1"]


##@memoized
def test_adem_2():
    assert adem_2(1,1) == {}
    assert adem_2(1,2) == {(3,) : 1}
    assert adem_2(2,2) == {(3,1) : 1}
    assert adem_2(4,2) == {(4, 2): 1} # Admissible
    assert adem_2(4,4) == {(6, 2): 1, (7, 1): 1}
    assert adem_2(5,7) == {}
    assert adem_2(6,7) == {(13,): 1, (12, 1): 1, (10, 3): 1}
    assert len(adem_2(100,100)) == 16
    assert len(adem_2(200,200)) == 23

#@memoized
def test_adem_generic():
    assert adem_generic(1,0,1, p = 3 ) == {(0, 2, 0) : 2}
    assert adem_generic(1,1,1, p = 3 ) == {(1, 2, 0) : 1, (0,2,1) : 1}
    assert adem_generic(1,0,2, p = 3 ) == {}
    assert adem_generic(1,1,2, p = 3 ) == {(0, 3, 1): 1, (1, 3, 0): 2}
    assert adem_generic(2,0,1, p = 3 ) == {}
    assert adem_generic(3,0,1, p = 3 ) == {(0, 3, 0, 1, 0): 1} # Admissible
    assert adem_generic(5,0,7, p = 3 ) == {(0, 11, 0, 1, 0): 1}
    assert adem_generic(25,1,20, p = 3) ==  {(1, 38, 0, 7, 0): 1, (0, 38, 1, 7, 0): 1, (0, 37, 1, 8, 0): 1}
    assert len(adem_generic(1000,1,1000, p = 3)) == 107
    
    assert adem_generic(1, 1, 2, p=5) == {(0, 3, 1): 1, (1, 3, 0): 2}
    
    assert adem_generic(1,0,1, p = 23) == {(0, 2, 0) : 2}    
    assert adem_generic(1,1,1, p = 23) == {(1, 2, 0) : 1, (0,2,1) : 1}    
    assert adem_generic(2,0,1, p = 23) == {(0,3,0) : 3}
    assert adem_generic(5,0,7, p = 23) == {(0, 12, 0): 10}
    assert len(adem_generic(2000,1,5000, p = 23)) == 9

def test_adem():
    """
        This is a dispatch to adem_generic or adem_2 so no need to test I don't think.
    """
    pass

#@memoized
def test_make_mono_admissible_2():
    """
        mono
        TODO: fill me in!
        SteenrodAlgebra(p=2, basis='adem').Q(2) * (Sq(6) * Sq(2)) # indirect doctest
        Sq^10 Sq^4 Sq^1 + Sq^10 Sq^5 + Sq^12 Sq^3 + Sq^13 Sq^2
    """
    assert make_mono_admissible_2((12,)) == {(12,): 1} # already admissible
    assert make_mono_admissible_2((2,1)) == {(2, 1): 1} # already admissible
    assert make_mono_admissible_2((2,2)) == {(3, 1): 1}
    assert make_mono_admissible_2((2, 2, 2)) == {(5, 1): 1}
    
#@memoized
def test_make_mono_admissible_generic():
    """
        mono, p
        TODO: add more test cases!
    """
    assert make_mono_admissible_generic((0, 2, 0, 1, 0), p=7) == {(0, 3, 0): 3}
    


def test_make_mono_admissible():
    """
        A dispatch to make_mono_admissible_2 or make_mono_admissible_generic.
    """
    pass

def test_product_2():
    """
        Just calls make_mono_admissible_2 on the concatenation of the args.
    """
    pass

def test_product_generic():
    """
        Probably should have a couple tests here to make sure concatenation is handled correctly.
    """
    pass

def test_product():
    """
        A dispatch to product_2 or product_generic.
    """
    pass


#@memoized
def test_basis_2():
    assert set(basis_2(0)) == set([()])
    assert set(basis_2(1)) == set([(1,)])
    assert set(basis_2(2)) == set([(2,)])
    assert set(basis_2(3)) == set([(3,),(2,1)])
    assert set(basis_2(4)) == set([(4,),(3,1)])
    assert set(basis_2(7)) == set([(4,2,1),(5,2), (6,1),(7,)])
    assert len(basis_2(100)) == 1189
    
#@memoized
def test_basis_generic():
    assert set(basis_generic(0, p = 3)) == set([()])
    assert set(basis_generic(1, p = 3)) == set([(1,)])
    assert set(basis_generic(2, p = 3)) == set([])
    assert set(basis_generic(4, p = 3)) == set([(0,1,0)])
    assert set(basis_generic(5, p = 3)) == set([(1,1,0),(0,1,1)])
    assert set(basis_generic(10, p = 3)) == set([(1,2,1)])
    assert len(basis_generic(400, p = 3)) == 395


def test_basis():
    """
        Another dispatch...
    """
    pass
