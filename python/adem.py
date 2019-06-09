"""
    File: adem.py
    Author: Hood Chatham

    Defines the basic operations on Adem algebras. No object oriented code,
    all functions take basis vectors as inputs and output dictionaries representing
    Fp-linear combinations of basis vectors. The relevant calls will be wrapped
    by methods of the AdemAlgebra and AdemElement classes in steenrod.py.
    They are systemtically extended over vector inputs using the
    @linearextension decorator.

    The algorithms here are all copied from Sage, in particular from the two
    files steenrod_algebra_mult.py and (the basis code from) steenrod_algebra_basis.py
    I have made significant improvements in code legibility and removed all references to Sage code.
"""


import math

import combinatorics
from memoized import memoized

class MinimalAdemAlgebra:
    """A record object with fields p and generic
       Many of the methods take an "algebra" argument which needs to be filled by
       an object with p and generic defined. The calling code in steenrod.py
       has an AdemAlgebra instance that it hands in. This little record
       is for testing purposes.
    """
    def __init__(self, p, generic=None):
        self.p = p
        if generic is None:
            generic = p != 2
        self.generic = generic


def adem_basis_elt_generic_map(*, P_fn, b, basis_elt):
    """Map (P_fn, b) over basis_elt.
       Generic adem basis elements alternate [b_0, P_1, b_1, ..., P_n, b_n].
       Apply P_fn to P_i and if b_i is a 1 replace it with b, otherwise drop it.
    """
    Ps = [P_fn(P) for P in basis_elt[1::2]]
    bs = [b if epsilon else None for epsilon in basis_elt[ ::2]]
    result = [None] * len(basis_elt)
    result[1::2] = Ps
    result[ ::2] = bs
    return [x for x in result if x is not None]


def adem_basis_elt_2_map(*, Sq_fn, basis_elt):
    """Map Sq_fn over basis_elt
       We don't actually use this since it seems silly.
    """
    return [Sq_fn(Sq) for Sq in basis_elt]


@memoized
def adem_2(a, b):
    """Return the adem relation Sqa * Sqb when p=2"""
    if b == 0:
        return {(a,) : 1}
    if a == 0:
        return {(b,) : 1}
    if a >= 2*b:
        return {(a, b) : 1}
    result = {}
    for j in range(1 + a//2):
        if combinatorics.binomial_2(b-j-1, a-2*j) == 1:
            if j == 0:
                result[(a+b,)] = 1
            else:
                result[(a+b-j, j)] = 1
    return result

@memoized
def adem_generic(A, bockstein, B, *, p):
    """Return the generic adem relation for P(A)*P(B) or P(A) * beta * P(B)"""
    if A == 0:
        return {(bockstein, B, 0) : 1}
    
    if B == 0:
        return {(0, A, bockstein) : 1}
    
    if A >= p*B + bockstein: # admissible
        return {(0, A, bockstein, B, 0) : 1}
    
    result = {}
    for e1 in range(1 + bockstein):
        e2 = bockstein - e1
        for j in range(1 + A//p):
            coeff = combinatorics.binomial_odd((B-j) * (p-1) - 1 + e1, A - p*j - e2, p)
            coeff *= (-1)**(A+j + e2)
            coeff = coeff % p
            if coeff != 0 and j == 0:
                result[(e1, A+B, e2)] = coeff
            elif coeff != 0 and j != 0:
                result[(e1, A+B-j, e2, j, 0)] = coeff 
    return result

def adem(a, b, c=None, *, algebra):
    r"""Get the Adem relation.

    This is a dispatch to adem_2 and adem_generic. It's never used in the code
    because all of the calling code also breaks into the generic case or not.

    Here's stolen Sage documentation:
    INPUT:

    - `a`, `b`, `c` (optional) - nonnegative integers, corresponding
      to either `P^a P^b` or (if `c` present) to `P^a \beta^b P^c`
    - algebra -- the algebra we're in

    OUTPUT:

    a dictionary representing the mod p Adem relation
    satisfied by P^a P^b or (if c present) P^a \beta^b P^c.

    If P^a P^b or P^a \beta^b P^c is admissible, then it is returned as is.
    Otherwise the right hand side of the Adem relation is returned.

    If the algebra is not generic, then c must be zero.

    If the algebra is generic, then with two inputs `a` and `b`, the
    Adem relation for `P^a P^b` is computed.  With three inputs `a`,
    `b`, `c`, the Adem relation for `P^a \beta^b P^c` is computed.
    In either case, the keys in the output are all tuples of odd length,
    with ``(i_1, i_2, ..., i_m)`` representing

        \beta^{i_1} P^{i_2} \beta^{i_3} P^{i_4} ... \beta^{i_m}
    """
    if algebra.generic:
        return adem_generic(a, b, c, p=algebra.p)
    else:
        if c is not None:
            raise ValueError("When p = 2, c should be None")
        return adem_2(a, b)


def make_mono_admissible_2(mono):
    result = {}
    make_mono_admissible_2_helper(result, list(mono))
    return result

def make_mono_admissible_2_helper(result, mono):
    """Reduce a monomial into a linear combination of admissible monomials when p = 2"""
    # check to see if admissible:
    nonadmissible_indices = [j for j in range(len(mono) - 1) if mono[j] < 2*mono[j+1]]
    if not nonadmissible_indices:
        mono = tuple(mono)
        if mono in result:
            del result[mono]
        else:
            result[mono] = 1
        return
    j = nonadmissible_indices[0]
    ab = mono[j:j+2]
    y = adem_2(mono[j], mono[j+1])
    for x in y:
        if len(x) == 2:
            mono[j:j+2] = x # modify mono in place =)
            new = mono
        else:
            new = mono[:j] + [x[0]] + mono[j+2:] # make a new shorter one.
        make_mono_admissible_2_helper(result, new)
    mono[j:j+2] = ab # Put monomial back the way we found it.

def make_mono_admissible_generic(mono, p):
    result = {}
    make_mono_admissible_generic_helper(result, 1, list(mono), p)
    return result

def make_mono_admissible_generic_helper(result, coeff, mono, p):
    """Reduce a monomial into a linear combination of admissible monomials for the generic Steenrod algebra"""
    # check to see if admissible:
    # print(mono, coeff)
    nonadmissible_indices = [
        j for j in range(1, len(mono) - 2, 2)
        if mono[j] < mono[j+1] + p * mono[j+2]
    ]
    if not nonadmissible_indices:
        # It's admissible, just add it to the result.
        mono = tuple(mono)
        if mono in result:
            result[mono] += coeff # y[x] * new[m]
            result[mono] = result[mono] % p
            if result[mono] == 0:
                del result[mono]
        else:
            result[mono] = coeff # * new[m]        
        return 

    j = nonadmissible_indices[0]
    y = adem_generic(*mono[j:j+3], p=p)
    for x in y:
        new_x = list(x)
        new_x[0] = mono[j-1] + x[0]
        # add bocksteins together
        if len(mono) >= j+3:
            new_x[-1] = mono[j+3] + x[-1]
        # If this next conditional fails, there are two bocksteins in a row.
        if new_x[0] <= 1 and new_x[-1] <= 1:
            new = mono[:j-1] + new_x + mono[j+4:]
            make_mono_admissible_generic_helper(result, y[x]*coeff, new, p)


def make_mono_admissible(mono, *, algebra):
    r"""Reduce a monomial into a linear combination of admissible monomials.
    This is a dispatch to make_mono_admissible_generic or make_mono_admissible_2.

    Stolen Sage documentation:

    Given a tuple ``mono``, view it as a product of Steenrod
    operations, and return a dictionary giving data equivalent to
    writing that product as a linear combination of admissible
    monomials.

    When p = 2 and the algebra is not generic, the sequence
    `(i_1, i_2, ...)` is admissible if `i_j \geq 2 i_{j+1}` for all `j`.

    When `p` is odd, the sequence `(e_1, i_1, e_2, i_2, ...)` is
    admissible if `i_j \geq e_{j+1} + p i_{j+1}` for all `j`.

    INPUT:

    - ``mono`` - a tuple of non-negative integers
    - algebra  - a class with properies p and generic

    OUTPUT:

    A linear combination of admissible monomials represented by a dictionary of
    terms of the form tuple: coeff, where 'tuple' is an admissible tuple of non-negative integers and
    'coeff' is its coefficient.    When `p` is odd, each tuple
    must have an odd length: it should be of the form `(e_1, i_1, e_2,
    i_2, ..., e_k)` where each `e_j` is either 0 or 1 and each `i_j`
    is a positive integer: this corresponds to the admissible monomial

    .. MATH::

       \beta^{e_1} \mathcal{P}^{i_2} \beta^{e_2} \mathcal{P}^{i_2} ...
       \mathcal{P}^{i_k} \beta^{e_k}

    ALGORITHM:

    Given `(i_1, i_2, i_3, ...)`, apply the Adem relations to the first
    pair (or triple when `p` is odd) where the sequence is inadmissible,
    and then apply this function recursively to each of the resulting
    tuples `(i_1, ..., i_{j-1}, NEW, i_{j+2}, ...)`, keeping track of
    the coefficients.
    """
    if len(mono) == 1:
        return {mono: 1}
    if not algebra.generic and len(mono) == 2:
        return adem_2(*mono)
    if algebra.generic:
        return make_mono_admissible_generic(mono, algebra.p)
    else:
        return make_mono_admissible_2(mono)

def product_2(m1, m2):
    """Multiply monomials m1 and m2 and write the result in the Adem basis for p = 2."""
    return make_mono_admissible_2(list(m1) + list(m2))    

def product_generic(m1, m2, p):
    """Multiply monomials m1 and m2 and write the result in the Adem basis in the generic case."""
    if m1[-1] == m2[0] == 1:
        return {}
    else:
        return make_mono_admissible_generic(m1[:-1] + (m1[-1] + m2[0],) + m2[1:], p)

def product(m1, m2, *, algebra):
    """Multiply monomials m1 and m2 and write the result in the Adem basis.
       This is a dispatch to product_2 and product_generic.
    """
    if algebra.generic:
        return product_generic(m1, m2, algebra.p)
    else:
        return product_2(m1, m2)




@memoized
def basis_2(n, *, bound=1):
    """Get the basis for the n dimensional part of the Steenrod algebra.
       Build basis recursively.  last = last term.
       last is >= bound, and we will append (last,) to the end of
       elements from basis (n - last, bound=2 * last).
       This means that 2 last <= n - last, or 3 last <= n.
    """
    if n == 0:
        return ((),)
    result = [(n,)]
    for last in range(bound, 1 + n // 3):
        for vec in basis_2(n - last, bound=2 * last):
            result.append(vec + (last,))
    return tuple(result)

@memoized
def basis_generic(n, *, p, bound=1):
    """Get the basis for the n dimensional part of the Steenrod algebra."""
    if n == 0:
        return ((0,),) # 
    if n == 1:
        return ((1,),)
    result = []
        
    # append P^{last} beta^{epsilon}
    for epsilon in [0,1]:
        # Without this lower bound edge case we lose the element (0, 1, 1) in degree 5.
        # I don't have a good explanation for what it means yet.
        lower_bound = bound + epsilon if bound > 1 else 1
        for last in range(lower_bound, 1 + (n // (2*(p - 1)))):
            remaining_degree = n - 2*(p-1)*last - epsilon
            basis_in_remaining_degree = basis_generic(remaining_degree, p=p, bound=p * last)
            for vec in basis_in_remaining_degree:
                result.append(vec + (last, epsilon))
    return tuple(result)

def basis(n, *, algebra):
    """Adem basis in dimension `n`. Dispatch to basis_2 and basis_generic.

    INPUT:
    - ``n`` - non-negative integer
    - ``algebra`` - an object with properties p and generic

    OUTPUT: tuple of mod p Serre-Cartan basis elements in dimension n
    """
    if algebra.generic:
        return basis_generic(n, p=algebra.p)
    else:
        return basis_2(n)


if __name__ == "__main__":
    dim = 0
    for x in basis_2(500):
        dim+=1
    print(dim)
