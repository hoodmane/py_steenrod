import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,os.pardir,"python"))

from milnor import *

#class FullProfile:

#class Profile:

def test_profile():
    # TODO: fill me in
    pass  

#class MinimalMilnorAlgebra:
   
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
    assert product_even((1,),(1,), 3) == {(2,) : 2}
    assert product_even((1,),(0, 1,), 3) == {(1,1) : 1}
    assert product_even((0,2), (1,),2) == {(1, 2): 1, (0, 0, 1): 1}
    assert product_even((0,3), (1,),3) == {(1, 3): 1, (0, 0, 1): 1}
    assert product_even((0,9), (1,1),3) == {(1, 10): 1, (0, 7, 1): 1, (1, 0, 0, 1): 1}

def test_product_full_Qpart():
    assert product_full_Qpart(((),(1,)),(0,), 3) == {((0,),(1,)) : 1, ((1,),()): 1}
    

def test_product_full():
    assert product_full(((),(1,)),((0,),()), 3) == {((0,),(1,)) : 1, ((1,),()): 1}
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
    assert set(basis_even(2, 2, Profile())) == set([(2,)])
    assert set(basis_even(3, 2, Profile())) == set([(0, 1), (3,)])
    assert set(basis_even(4, 2, Profile())) == set([(1, 1), (4,)])
    assert set(basis_even(4, 2, Profile(profile=[2,1]))) == set([(1, 1)])
    assert set(basis_even(4, 2, Profile(profile=(), truncation=0))) == set([])
    assert set(basis_even(4, 2, Profile(profile=(), truncation=Infinity))) == set([(1, 1), (4,)])
    assert set(basis_even(7, 2, Profile())) == set([(0, 0, 1), (1, 2), (4, 1), (7,)])


def test_basis_generic_Q_part():
    """
        q_deg, p, profile
    """
    pass

def test_basis_generic():
    """
        n, p, profile
    """
    assert set(basis_generic(1, 3, FullProfile())) == set([((0,), ())])
    assert set(basis_generic(9, 3, FullProfile())) == set([((1,), (1,)), ((0,), (2,))])
    assert set(basis_generic(17, 3, FullProfile())) == set([((2,), ()), ((1,), (3,)), ((0,), (0, 1)), ((0,), (4,))])
    assert set(basis_generic(48, 5, FullProfile())) == set([((), (0, 1)), ((), (6,))])
    assert len(basis_generic(100, 3, FullProfile())) == 13
    assert len(basis_generic(200, 7, FullProfile())) == 0
    assert len(basis_generic(240, 7, FullProfile())) == 3
    assert len(basis_generic(240, 7, FullProfile(odd_part = (), even_part = (), truncation = Infinity))) == 3
    assert len(basis_generic(240, 7, FullProfile(odd_part = (), even_part = (), truncation = 0))) == 0    

def test_basis():
    """
        Just a dispatch
    """
    pass
        


