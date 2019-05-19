"""
    File: milnor.py
    Author: Hood Chatham

    Defines the basic operations on Milnor algebras. No object oriented code,
    all functions take basis vectors as inputs and output dictionaries representing
    Fp-linear combinations of basis vectors. The relevant calls will be wrapped
    by methods of the MilnorAlgebra and MilnorElement classes in steenrod.py.
    They are systemtically extended over vector inputs using the
    @linearextension decorator.

    The algorithms here are all copied from Sage, in particular from the two
    files steenrod_algebra_mult.py and (the basis code from) steenrod_algebra_basis.py
    I have made significant improvements in code legibility and removed all references to Sage code.
"""

import itertools

import combinatorics
from memoized import memoized
from infinity import Infinity

class FullProfile:
    """A "Full" profile is the generic case. It has an "even" Profile and an "odd" Profile.
       The "even" profile restricts the P's and the "odd" part restricts the Q's.
    """
    def __init__(self, odd_part=None, even_part=None, truncation=None):
        self.even_restricted = not (even_part is None and truncation is None)
        self.odd_restricted = not (odd_part is None and truncation is None)
        self.restricted = self.even_restricted or self.odd_restricted
        self.even_part = even_part or []
        self.odd_part = odd_part or []
        self.truncation = truncation
        if truncation is None and (even_part or odd_part):
            self.truncation = 0
        elif truncation is None:
            self.truncation = Infinity

    def get_even_part(self):
        return Profile(self.even_part, self.truncation)

    def get_odd_part(self):
        return Profile(self.odd_part, self.truncation)

class Profile:
    """ Encodes a Profile function.
        profile[i] returns the ith element of the profile list if present,
        otherwise it returns the truncation. This little class allows us to remove
        a large amount of conditional nonsense when checking whether a basis element
        is in the subalgebra.
    """
    def __init__(self, profile=None, truncation=None):
        self.even_restricted = not (profile is None and truncation is None)
        self.restricted = self.even_restricted
        self.profile = profile or []
        self.truncation = truncation
        if truncation is None and profile:
            self.truncation = 0
        elif truncation is None:
            self.truncation = Infinity

    def __getitem__(self, index):
        if callable(self.profile):
            return self.profile(index)
        if index < len(self.profile):
            return self.profile[index]
        else:
            return self.truncation

    def exponent(self, i, p):
        """Get p ^ profile[i].
           The reason we have this function is that we would prefer not to calculate
           p ^ infinity.
        """
        res = self[i]
        if res < Infinity:
            res = p ** res
        return res

class MinimalMilnorAlgebra:
    def __init__(self, p, generic=None, profile=None, truncation_type=None):
        self.p = p
        if generic is None:
            generic = p != 2
        self.generic = generic
        if generic:
            self.profile = FullProfile(
                profile and profile[0],
                profile and profile[1],
                truncation_type
            )
        else:
            self.profile = Profile(profile, truncation_type)

def initialize_milnor_matrix(r, s):
    """Initializes an len(r)+1 by len(s)+1 matrix
       Puts r along the first column and s along the first row and zeroes everywhere else.
    """
    rows = len(r) + 1
    cols = len(s) + 1
    M = [[0 for _ in range(cols)] for _ in range(rows)]
    for j in range(1, cols):
        M[0][j] = s[j-1]
    for i in range(1, rows):
        M[i][0] = r[i-1]
    return M

def step_milnor_matrix(M, r, s, i, j, x):
    """
        This seems to move an i x j block of M back to the first row and column.
        To be honest, I don't really know what the point is, but the milnor_matrices
        function was a little long and this seemed like a decent chunk to extract.
        At least it contains all of the steps that modify M so that seems like a good thing.
    """
    cols = len(s) + 1
    for row in range(1, i):
        M[row][0] = r[row-1]
        for col in range(1, cols):
            M[0][col] += M[row][col]
            M[row][col] = 0
    for col in range(1, j):
        M[0][col] += M[i][col]
        M[i][col] = 0
    M[0][j] -= 1
    M[i][j] += 1
    M[i][0] = x

def milnor_matrices(r, s, p):
    """
        Generator for Milnor matrices. milnor_product_even iterates over this.
        Uses the same algorithm Monks does in his Maple package to iterate through
        the possible matrices: see
        https://monks.scranton.edu/files/software/Steenrod/steen.html
    """
    rows = len(r) + 1
    cols = len(s) + 1
    M = initialize_milnor_matrix(r, s)
    yield M
    found = True
    while found:
        found = False
        for i in range(1, rows):
            if found:
                break
            total = M[i][0]
            for j in range(1, cols):
                p_to_the_j = p ** j
                # check to see if column index j is small enough or
                # there's nothing above the column to add to it
                if total < p_to_the_j or all([M[k][j] == 0 for k in range(i)]):
                    total += M[i][j] * p_to_the_j
                else:
                    step_milnor_matrix(M, r, s, i, j, total - p_to_the_j)
                    found = True
                    yield M
                    break

def remove_trailing_zeroes(l):
    """Remove trailing zeroes from the list l."""
    for i in range(len(l) - 1, -1, -1):
        if l[i] != 0:
            return None
        else:
            l.pop()

@memoized
def product_even(r, s, p):
    """Handles the multiplication in the even subalgebra of the Steenrod algebra P.
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
            nth_diagonal = [M[i][n-i] for i in range(max(0, n-cols+1), min(1+n, rows))]
            coeff *= combinatorics.multinomial(nth_diagonal, p)
            coeff = coeff % p
            if coeff == 0:
                break
            diagonal_sums[n-1] = sum(nth_diagonal)
        if coeff != 0:
            remove_trailing_zeroes(diagonal_sums)
            t = tuple(diagonal_sums)
            old_coeff = result[t] if (t in result) else 0
            result[t] = (coeff + old_coeff) % p
    return result

#@memoized
def product_full_Qpart(m1, f, p):
    """Reduce m1 * f = (Q_e0 Q_e1 ... P(r1, r2, ...)) * (Q_f0 Q_f1 ...) into the form Sum of Q's * P's
       Result is represented as dictionary of pairs of tuples.
    """
    result = {m1: 1}
    for k in f:
        old_result = result
        print(old_result)
        result = {}
        p_to_the_k = p**k
        for mono in old_result: 
            for i in range(0, 1 + len(mono[1])):
                q_mono = mono[0]
                p_mono = mono[1]
                if (k + i in q_mono):
                    continue
                if (i > 0 and p_mono[i - 1] < p_to_the_k):
                    continue

                if i > 0:
                    p_mono = list(p_mono)
                    p_mono[i - 1] -= p_to_the_k
                    remove_trailing_zeroes(p_mono)
                    p_mono = tuple(p_mono)                

                ind = len([x for x in q_mono if x >= k+i])                
                q_mono = q_mono[:len(q_mono) - ind] + (k+i,) + q_mono[len(q_mono) - ind:]
                
                coeff = (-1)**ind * old_result[mono]
                result[(q_mono, p_mono)] = coeff % p
    return result

#@memoized
def product_full(m1, m2, p):
    r"""
    Product of Milnor basis elements defined by m1 and m2 at the prime p.

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
        sorted(milnor.product_generic(((0,2),(5,)), ((1,),(1,)), 5).items())
        [(((0, 1, 2), (0, 1)), 4), (((0, 1, 2), (6,)), 4)]
        milnor.product_generic(((0,2,4),()), ((1,3),()), 7)
        {((0, 1, 2, 3, 4), ()): 6}
        milnor.product_generic(((0,2,4),()), ((1,5),()), 7)
        {((0, 1, 2, 4, 5), ()): 1}
        sorted(milnor.product_generic(((),(6,)), ((),(2,)), 3).items())
        [(((), (0, 2)), 1), (((), (4, 1)), 1), (((), (8,)), 1)]

        Associativity once failed because of a sign error:

        sage: a,b,c = A.Q_exp(0,1), A.P(3), A.Q_exp(1,1)
        sage: (a*b)*c == a*(b*c)
        True
    """
    (f, s) = m2
    m1_times_f = product_full_Qpart(m1, f, p)
    # Now for the Milnor matrices.  For each entry '(e,r): coeff' in answer,
    # multiply r with s.  Record coefficient for matrix and multiply by coeff.
    # Store in 'result'.
    if not s:
        return m1_times_f
    result = {}
    for (e, r) in m1_times_f:
        coeff = m1_times_f[(e, r)]
        prod = product_even(r, s, p=p)
        for k in prod:
            old_coeff = result[(e, k)] if (e, k) in result else 0
            result[(e, k)] = (old_coeff + coeff*prod[k]) % p
    return result

def product_2(r, s):
    """Multiplication of Milnor basis elements in the non generic case."""
    return product_even(r, s, 2)

product_generic = product_full

def product(r, s, *, algebra):
    """Multiply r and s in the Milnor algebra determined by algebra.
       Note that since profile functions determine subalgebras, the product
       doesn't need to care about the profile function at all.
    """
    if algebra.generic: 
        return product_full(r, s, algebra.p)
    else:
        return product_even(r, s, algebra.p)


def check_even_profile(p, profile, exponents):
    if not profile.even_restricted:
        return True
    for (i, exp) in enumerate(exponents):
        if exp >= profile.exponent(i, p):
            return False
    return True
    
def check_odd_profile(profile, q_mono):
    if not profile.restricted:
        return True
    for i in q_mono:
        if profile[i] <= 1:
            return False
    return True

def basis_even(n, p, profile):
    """Return the even part of the basis in degree n * 2*(p-1).
       In the nongeneric case, this actually just gets the whole degree n basis.
       Note the factor of two difference between 2*(2-1) and 1.
    """
    if n == 0:
        return [()]
    result = []
    for mono in combinatorics.WeightedIntegerVectors(n, combinatorics.xi_degrees(n, p=p)):
        exponents = list(mono)
        remove_trailing_zeroes(exponents)
        if check_even_profile(p, profile, exponents):
            result.append(tuple(exponents))
    return tuple(result)

def basis_generic_Q_part(q_deg, p, profile):
    """Returns the "Q-part" of the basis in degree q_deg.
       This means return the set of monomials Q(i_1) * ... * Q(i_k) where i_1 < ... < i_k
       and the product is in q_deg. Basically it's just an issue of finding partitions of
       q_deg into parts of size |Q(i_j)|, and then there's a profile condition.
    """
    q_degrees = combinatorics.tau_degrees(q_deg, p=p)
    result = []
    for sigma in combinatorics.restricted_partitions(q_deg, q_degrees):
        # q_mono is the list of indices ocurring in the partition
        q_mono = [idx for idx in range(len(sigma)) if sigma[idx] == 1]
        if check_odd_profile(profile, q_mono):
            result.append(tuple(q_mono))
    return result

def basis_generic(n, p, profile):
    """Get the basis in degree n for the generic steenrod algebra at the prime p.
       Basically we just put together the "even part" of the basis and the "Q part".
    """
    if n == 0:
        return (((), ()),)
    result = []
    # p_deg records the desired degree of the P part of the basis element.
    # Since p-parts are always divisible by 2p-2, we divide by this first.
    q = 2*(p-1)
    min_q_deg = (-1 + p ** (n % q - 1)) // (p - 1)
    for p_deg in range(n // (2 * (p - 1)) + 1):
        q_deg = n - 2 * p_deg * (p - 1)

        # if this inequality holds, no way to have a partition
        # with distinct parts.
        if q_deg < min_q_deg:
            break

        Q_parts = basis_generic_Q_part(q_deg, p, profile.get_odd_part())
        P_parts = basis_even(p_deg, p, profile.get_even_part())
        for q_p_pair in itertools.product(Q_parts, P_parts):
            result.append(q_p_pair)
    return tuple(result)


def basis(n, *, algebra):
    r"""
    Milnor basis in dimension `n` with profile function ``profile``.

    INPUT:

    - ``n`` - non-negative integer

    - ``p`` - positive prime number

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
    if algebra.generic:
        return basis_generic(n, algebra.p, algebra.profile)
    else:
        return basis_even(n, 2, algebra.profile)
