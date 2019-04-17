from itertools import product

import combinatorics
from memoized import memoized
from infinity import Infinity

def unit(p, generic = None):
    if p > 2 or generic :
        return { ((),()) : 1 }
    else:
        return { () : 1 }
        
        
def ademSq(i):
    return { (i,) : 1 } 
    
def ademP(i):
    return { ((i,),()) : 1 }
    
def Q(i):
    return { ((),(i,)) : 1 }
    

def step_milnor_matrix(M, i, j):
    """
        This seems to move an i x j block of M back to the first row and column.
        To be honest, I don't really know what the point is, but the milnor_matrices 
        function was a little long and this seemed like a decent chunk to extract.
    """
    for row in range(1,i):
        M[row][0] = r[row-1]
        for col in range(1,cols):
            M[0][col] += M[row][col]
            M[row][col] = 0
    for col in range(1,j):
        M[0][col] += M[i][col]
        M[i][col] = 0
    M[0][j] -= 1
    M[i][j] += 1 

def milnor_matrices(r, s, p):
    r"""
        Generator for Milnor matrices. milnor_product_even iterates over this.
        Uses the same algorithm Monks does in his Maple package to iterate through the possible matrices: see
        https://monks.scranton.edu/files/software/Steenrod/steen.html
    """
    rows = len(r) + 1;
    cols = len(s) + 1;
    diags = len(r) + len(s)
    # initialize matrix
    M = [[0 for col in range(cols)] for row in range(rows)]
    for j in range(1,cols):
        M[0][j] = s[j-1]
    for i in range(1,rows):
        M[i][0] = r[i-1]
    yield M
    found = True
    while found:
        found = False
        for i in range(1,rows):
            if found:
                break
            sum = M[i][0]
            for j in range(1, cols):
                p_to_the_j = p ** j
                # check to see if column index j is small enough or there's nothing above the column to add to it
                if sum < p_to_the_j or all([M[k][j] == 0 for k in range(i)]):
                    sum += M[i][j] * p_to_the_j
                else:
                    milnor_matrix_next(M, i, j)
                    M[i][0] = sum - p_to_the_j 
                    found = True
                    yield M
                    break

@memoized
def product_even(r, s, p):
    r"""
        Handles the multiplication in the even subalgebra of the Steenrod algebra P.
        When p = 2, this is isomorphic to the whole Steenrod algebra so this method does everything.
    """
    result = {}
    rows = len(r) + 1
    cols = len(s) + 1
    diags = len(r) + len(s)
    for M in milnor_matrices(r, s, p):
        # check diagonals
        coeff = 1
        diagonal_sums = [0] * diags
        for n in range(1, diags + 1):
            nth_diagonal = [M[i][n-i] for i in range(max(0,n-cols+1), min(1+n,rows))]
            coeff *= combinatorics.multinomial(nth_diagonal, p)
            coeff = coeff % p
            if coeff == 0:
                break
            diagonal_sums[n-1] = sum(nth_diagonal)
        if coeff != 0:
            i = diags - 1
            while i >= 0 and diagonal_sums[i] == 0:
                i -= 1
            t = tuple(diagonal_sums[:i+1])
            old_coeff = result[t] if (t in result) else 0
            result[t] = (coeff + old_coeff) % p
    return result


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
                    if len(q_mono) > 0:
                        ind = len(filter(lambda x : x >= k+i, q_mono))
                    else:
                        ind = 0
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

@memoized
def product_full(m1,m2,p):
    r"""
    Product of Milnor basis elements defined by m1 and m2 at the odd prime p.

    INPUT:

    - m1 - pair of tuples (e,r), where e is an increasing tuple of
      non-negative integers and r is a tuple of non-negative integers
    - m2 - pair of tuples (f,s), same format as m1
    - p -- odd prime number

    OUTPUT:

    Dictionary of terms of the form (tuple: coeff), where 'tuple' is
    a pair of tuples, as for r and s, and 'coeff' is an integer mod p.

    This computes the product of the Milnor basis elements
    $Q_{e_1} Q_{e_2} ... P(r_1, r_2, ...)$ and
    $Q_{f_1} Q_{f_2} ... P(s_1, s_2, ...)$.

    EXAMPLES::
        import milnor
        sorted(milnor.product_odd(((0,2),(5,)), ((1,),(1,)), 5).items())
        [(((0, 1, 2), (0, 1)), 4), (((0, 1, 2), (6,)), 4)]
        milnor.product_odd(((0,2,4),()), ((1,3),()), 7)
        {((0, 1, 2, 3, 4), ()): 6}
        milnor.product_odd(((0,2,4),()), ((1,5),()), 7)
        {((0, 1, 2, 4, 5), ()): 1}
        sorted(milnor.product_odd(((),(6,)), ((),(2,)), 3).items())
        [(((), (0, 2)), 1), (((), (4, 1)), 1), (((), (8,)), 1)]    Associativity once failed because of a sign error::

        sage: a,b,c = A.Q_exp(0,1), A.P(3), A.Q_exp(1,1)
        sage: (a*b)*c == a*(b*c)
        True
    """
    (f,s) = m2
    m1_times_f = product_full_Qpart(m1, f, p)
    # Now for the Milnor matrices.  For each entry '(e,r): coeff' in answer,
    # multiply r with s.  Record coefficient for matrix and multiply by coeff.
    # Store in 'result'.
    if not s:
        result = m1_times_f
    else:
        result = {}
        for (e, r) in m1_times_f:
            coeff = m1_times_f[(e,r)]
            # Milnor multiplication for r and s
            prod = product_even(r, s, p)
            for k in prod:
                old_coeff = result[(e,k)] if (e, k) in result else 0
                result[(e,k)] = (old_coeff + coeff*prod[k]) % p
    return result

def product_2(r, s):
    return product_even(r, s, 2)

product_odd = product_full

def product(r, s, p, generic = None):
    if generic is None:
        generic = p != 2
    if generic: # Should also be not generic...
        return product_full(r, s, p)
    else:
        return product_even(r, s, p)
        
        
class FullProfile:
    def __init__(self, even_part = None, odd_part = None, truncation = None):
        self.even_part = even_part or []
        self.odd_part = odd_part or []
        self.truncation = truncation
        if truncation is None:
            if even_part or odd_part:
                self.truncation = 0
            else:
                self.truncation = Infinity
    
    def getEvenPart(self):
        return Profile(self.even_part, self.truncation)
    
    def getOddPart(self):
        return Profile(self.odd_part, self.truncation)

class Profile:
    def __init__(self, profile, truncation = None):
        self.profile = profile or []
        self.truncation = truncation
        if truncation is None:
            if profile:
                self.truncation = 0
            else:
                self.truncation = Infinity        

    def __getitem__(self,index):
        if callable(self.profile):
            return profile(index)
        if index < len(self.profile):
            return self.profile[index]
        else:
            return self.truncation
    
    def exponent(self, i, p):
        res = self[i]
        if res < Infinity:
            res = p ** res    
        return res        

def basis_even(n, p, profile):
    if n == 0:
        return [[]]
    result = []
    for mono in combinatorics.WeightedIntegerVectors(n, combinatorics.xi_degrees(n, p, reverse = False)):
        exponents = list(mono)
        while len(exponents) > 0 and exponents[-1] == 0: 
            exponents.pop()
        okay = True
        for i in range(len(exponents)):
            if exponents[i] >= profile.exponent(i, p):
                okay = False
                break
        if okay:
            result.append(tuple(exponents)) 
    return tuple(result)  

def basis_odd_Q_part(q_deg, p, profile):
    q_degrees = combinatorics.xi_degrees((q_deg - 1)//(2*(p-1)), p)
    q_degrees = [ 1+2*(p-1)*d for d in q_degrees ]
    q_degrees.append(1);
    q_degrees_decrease = list(q_degrees);
    q_degrees.reverse();
    result = []
    for sigma in combinatorics.restricted_partitions(q_deg, q_degrees_decrease, True):
        # q_mono is the list of indices ocurring in the partition
        q_mono = [idx for (idx, q_deg) in enumerate(q_degrees) if q_deg in sigma]
        # check profile:
        okay = True;
        for i in q_mono:
            if profile[i] <= 1:
                okay = False;
                break
        if okay:
            result.append(q_mono)
    return result

def basis_odd(n, p, profile):
    if n == 0:
        return (((),()))
    result = []
    # p_deg records the desired degree of the P part of the basis element.
    # Since p-parts are always divisible by 2p-2, we divide by this first.
    min_q_deg = (-1 + p ** ((n % (2 * (p - 1))) - 1)) // (p-1)
    for p_deg in range(n//(2*(p-1)) + 1):
        q_deg = n - 2 * p_deg * (p - 1);
        
        # if this inequality holds, no way to have a partition
        # with distinct parts.
        if q_deg < min_q_deg:
            break;
        
        P_parts = basis_even(p_deg, p, profile.getEvenPart())
        Q_parts = basis_odd_Q_part(q_deg, p, profile.getOddPart())
        for (p_mono, q_mono) in product(P_parts, Q_parts):
            result.append((tuple(q_mono), tuple(p_mono)));
    return tuple(result);
    

def basis(n, p=2, **kwds):
    r"""
    Milnor basis in dimension `n` with profile function ``profile``.

    INPUT:

    - ``n`` - non-negative integer

    - ``p`` - positive prime number (optional, default 2)

    - ``profile`` - profile function (optional, default None).
      Together with ``truncation_type``, specify the profile function
      to be used; None means the profile function for the entire
      Steenrod algebra.  See
      :mod:`sage.algebras.steenrod.steenrod_algebra` and
      :func:`SteenrodAlgebra <sage.algebras.steenrod.steenrod_algebra.SteenrodAlgebra>`
      for information on profile functions.

    - ``truncation_type`` - truncation type, either 0 or Infinity
      (optional, default Infinity if no profile function is specified,
      0 otherwise)

    OUTPUT: tuple of mod p Milnor basis elements in dimension n

    At the prime 2, the Milnor basis consists of symbols of the form
    `\text{Sq}(m_1, m_2, ..., m_t)`, where each
    `m_i` is a non-negative integer and if `t>1`, then
    `m_t \neq 0`. At odd primes, it consists of symbols of the
    form `Q_{e_1} Q_{e_2} ... P(m_1, m_2, ..., m_t)`,
    where `0 \leq e_1 < e_2 < ...`, each `m_i` is a
    non-negative integer, and if `t>1`, then
    `m_t \neq 0`.

    EXAMPLES::

        sage: import milnor
        sage: milnor.basis(7)
        ((0, 0, 1), (1, 2), (4, 1), (7,))
        sage: milnor.basis(7, 2)
        ((0, 0, 1), (1, 2), (4, 1), (7,))
        sage: milnor.basis(4, 2)
        ((1, 1), (4,))
        sage: milnor.basis(4, 2, profile=[2,1])
        ((1, 1),)
        sage: milnor.basis(4, 2, profile=(), truncation_type=0)
        ()
        sage: milnor.basis(4, 2, profile=(), truncation_type=Infinity)
        ((1, 1), (4,))
        sage: milnor.basis(9, 3)
        (((1,), (1,)), ((0,), (2,)))
        sage: milnor.basis(17, 3)
        (((2,), ()), ((1,), (3,)), ((0,), (0, 1)), ((0,), (4,)))
        sage: milnor.basis(48, p=5)
        (((), (0, 1)), ((), (6,)))
        sage: len(milnor.basis(100,3))
        13
        sage: len(milnor.basis(200,7))
        0
        sage: len(milnor.basis(240,7))
        3
        sage: len(milnor.basis(240,7, profile=((),()), truncation_type=Infinity))
        3
        sage: len(milnor.basis(240,7, profile=((),()), truncation_type=0))
        0
    """
    generic = kwds.get('generic', False if p==2 else True)
    from infinity import Infinity
    profile = kwds.get("profile", None)
    trunc = kwds.get("truncation_type", None)

    if not generic:
        profile = Profile(profile, trunc)
        return basis_even(n, 2, profile)
    else:  # p odd
        profile = FullProfile(profile and profile[0], profile and profile[1], trunc)
        return basis_odd(n, p, profile)
        


