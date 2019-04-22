import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,"src"))

from milnor import *

#class FullProfile:

#class Profile:

def test_profile():
    # TODO: fill me in
    pass  

class MinimalMilnorAlgebra:
    def __init__(self, p, generic = None, profile = None, truncation_type = None):
        self.p = p
        if generic is None:
            generic = p != 2
        self.generic = generic
        if generic:
            self.profile = FullProfile(profile and profile[0], profile and profile[1], truncation_type)
        else:  
            self.profile = Profile(profile, truncation_type)    
    
def test_initialize_milnor_matrix():
    M = initialize_milnor_matrix([4,5,6],["a","b","c", "d", "e"])
    assert len(M) == 4
    assert len(M[0]) == 6
    assert M[0] == [0, 'a', 'b', 'c', 'd', 'e']
    assert [r[0] for r in M] == [0,4,5,6]
    rest_of_M = [c[1:] for c in M[1:]]
    assert all([item==0 for sublist in rest_of_M for item in sublist])

def test_step_milnor_matrix():
    """
        This seems to move an i x j block of M back to the first row and column.
        To be honest, I don't really know what the point is, but the milnor_matrices 
        function was a little long and this seemed like a decent chunk to extract.
        
        I'm not really sure what the point of this method is, so I don't have any 
        test cases. I'll probably generate some later by running it on some inputs.
        M, r, s, i, j, x
    """
    pass

def test_milnor_matrices():
    r"""
        
    """
    # Since the generator writes each matrix into the same memory location, we
    # need to copy the matrices in order to list them.
    def mm_gen_to_list(generator):
        return [[c[:] for c in M] for M in generator]
    assert mm_gen_to_list(milnor_matrices([],[],2)) == [[[0]]]
    assert mm_gen_to_list(milnor_matrices([],[5],2)) == [[[0,5]]]
    assert mm_gen_to_list(milnor_matrices([5],[],2)) == [[[0],[5]]]
    assert mm_gen_to_list(milnor_matrices([1],[1],2)) == [[[0,1],[1,0]]]
    assert mm_gen_to_list(milnor_matrices([0,2],[1],2)) == [[[0, 1], [0, 0], [2, 0]], [[0, 0], [0, 0], [0, 1]]]
    assert mm_gen_to_list(milnor_matrices([0,4],[1,1],2)) == [[[0, 1, 1], [0, 0, 0], [4, 0, 0]], [[0, 0, 1], [0, 0, 0], [2, 1, 0]], [[0, 1, 0], [0, 0, 0], [0, 0, 1]]]
        
    assert mm_gen_to_list(milnor_matrices([0,3],[1],3)) == [[[0, 1], [0, 0], [3, 0]], [[0, 0], [0, 0], [0, 1]]]
    assert mm_gen_to_list(milnor_matrices([0,9],[1,1],3)) == [[[0, 1, 1], [0, 0, 0], [9, 0, 0]], [[0, 0, 1], [0, 0, 0], [6, 1, 0]], [[0, 1, 0], [0, 0, 0], [0, 0, 1]]]
        

def test_remove_trailing_zeroes():
    l = [1,0,0,1,0,0,1]; remove_trailing_zeroes(l)
    assert l == [1,0,0,1,0,0,1]
    l = [1,0,0,1,0,0,0]; remove_trailing_zeroes(l)
    assert l == [1,0,0,1]
    l = [1,0,0,0,0,0,0]; remove_trailing_zeroes(l)
    assert l == [1]
    l = [0,0,0,0,0,0,0]; remove_trailing_zeroes(l)
    assert l == []
    l = []; remove_trailing_zeroes(l)
    assert l == []

def test_product_even():
    r"""
        Handles the multiplication in the even subalgebra of the Steenrod algebra P.
        When p = 2, this is isomorphic to the whole Steenrod algebra so this method does everything.
    """
    


def product_full_Qpart(m1, f, p):
    """
        Reduce m1 * f = (Q_e0 Q_e1 ... P(r1, r2, ...)) * (Q_f0 Q_f1 ...) into the form Q's * P's
        Result is represented as dictionary of pairs of tuples.
    """
    result = {m1: 1}
    for k in f:
        old_result = result
        result = {}
        p_to_the_k = p**k
        for mono in old_result:
            for i in range(0,1+len(mono[1])):
                if (k+i not in mono[0]) and (i == 0 or p_to_the_k <= mono[1][i-1]):
                    q_mono = set(mono[0])
                    ind = len([ x for x in q_mono if x >= k+i ])
                    coeff = (-1)**ind * old_result[mono]
                    lst = list(mono[0])
                    if ind == 0:
                        lst.append(k+i)
                    else:
                        lst.insert(-ind,k+i)
                    q_mono = tuple(lst)
                    p_mono = list(mono[1])
                    if i > 0:
                        p_mono[i-1] = p_mono[i-1] - p_to_the_k
                    # The next two lines were added so that p_mono won't
                    # have trailing zeros. This makes p_mono uniquely
                    # determined by P(*p_mono).
                    while len(p_mono) > 0 and p_mono[-1] == 0:
                        p_mono.pop()
                    result[(q_mono, tuple(p_mono))] = coeff % p
    return result

def test_product_full():
    r"""
    
    """
    assert product_full(((0,2),(5,)), ((1,),(1,)), 5)== {((0, 1, 2), (0, 1)) :  4, ((0, 1, 2), (6,)): 4 }
    assert product_full(((0,2,4),()), ((1,3),()), 7) == {((0, 1, 2, 3, 4), ()): 6}
    assert product_full(((0,2,4),()), ((1,5),()), 7) == {((0, 1, 2, 4, 5), ()): 1}
    assert product_full(((),(6,)), ((),(2,)), 3) == {((), (0, 2)): 1, ((), (4, 1)): 1, ((), (8,)): 1}
    

def test_product_2():
    """
        TODO: test!
    """
    pass
            

def test_basis_even():
    """
        n, p, profile
    """
    assert set(basis_even(7, 2, Profile())) == set([(0, 0, 1), (1, 2), (4, 1), (7,)])
    assert set(basis_even(4, 2, Profile())) == set([(1, 1), (4,)])
    assert set(basis_even(4, 2, Profile(profile=[2,1]))) == set([(1, 1)])
    assert set(basis_even(4, 2, Profile(profile=(), truncation=0))) == set([])
    assert set(basis_even(4, 2, Profile(profile=(), truncation=Infinity))) == set([(1, 1), (4,)])

def test_basis_generic_Q_part():
    """
        q_deg, p, profile
    """
    pass

def test_basis_generic():
    """
        n, p, profile
    """
    assert set(basis_generic(9, 3, FullProfile())) == set([((1,), (1,)), ((0,), (2,))])
    assert set(basis_generic(17, 3, FullProfile())) == set([((2,), ()), ((1,), (3,)), ((0,), (0, 1)), ((0,), (4,))])
    assert set(basis_generic(48, 5, FullProfile())) == set([((), (0, 1)), ((), (6,))])
    assert len(basis_generic(100, 3, FullProfile())) == 13
    assert len(basis_generic(200, 7, FullProfile())) == 0
    assert len(basis_generic(240, 7, FullProfile())) == 3
    assert len(basis_generic(240, 7, FullProfile(even_part = (), odd_part = (), truncation = Infinity))) == 3
    assert len(basis_generic(240, 7, FullProfile(even_part = (), odd_part = (), truncation = 0))) == 0    

def test_basis():
    """
        Just a dispatch
    """
    pass
        


