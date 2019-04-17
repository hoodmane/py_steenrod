from __future__ import division

from math import ceil as ceiling

import combinatorics
from memoized import memoized

def unit():
    return { () : 1 }

@memoized
def adem_2(a, b, c=0):
    if b == 0:
        return {(a,): 1}
    elif a == 0:
        return {(b,): 1}
    elif a >= 2*b:
        return {(a,b): 1}
    result = {}
    for c in range(1 + a//2):
        if combinatorics.binomial_mod2(b-c-1, a-2*c) == 1:
            if c == 0:
                result[(a+b,)] = 1
            else:
                result[(a+b-c,c)] = 1
    return result

@memoized
def adem_odd(a, b, c, p):
    if a == 0 and b == 0:
            return {(c,): 1}
    if c == 0:
        bockstein = 0
        A = a
        B = b
    else:
        A = a
        B = c
        bockstein = b # should be 0 or 1
    if A == 0:
        return {(bockstein, B, 0): 1}
    if B == 0:
        return {(0, A, bockstein): 1}
    if bockstein == 0:
        if A >= p*B: # admissible
            return {(0,A,0,B,0): 1}
        result = {}
        for j in range(1 + a//p):
            coeff = (-1)**(A+j) * combinatorics.binomial_modp((B-j) * (p-1) - 1, A - p*j, p)
            if coeff % p != 0:
                if j == 0:
                    result[(0,A+B,0)] = coeff
                else:
                    result[(0,A+B-j,0,j,0)] = coeff
    else:
        if A >= p*B + 1: # admissible
            return {(0,A,1,B,0): 1}
        result = {}
        for j in range(1 + a//p):
            coeff = (-1)**(A+j) * combinatorics.binomial_modp((B-j) * (p-1), A - p*j, p)
            if coeff % p != 0:
                if j == 0:
                    result[(1,A+B,0)] = coeff
                else:
                    result[(1,A+B-j,0,j,0)] = coeff
        for j in range(1 + (a-1)//p):
            coeff = (-1)**(A+j-1) * combinatorics.binomial_modp((B-j) * (p-1) - 1, A - p*j - 1, p)
            if coeff % p != 0:
                if j == 0:
                    result[(0,A+B,1)] = coeff
                else:
                    result[(0,A+B-j,1,j,0)] = coeff
    return result


def adem(a, b, c=0, p=2, generic=None):
    r"""
    The mod `p` Adem relations

    INPUT:

    - `a`, `b`, `c` (optional) - nonnegative integers, corresponding
      to either `P^a P^b` or (if `c` present) to `P^a \beta^b P^c`
    - `p` - positive prime number (optional, default 2)
    - `generic` - whether to use the generic Steenrod algebra, (default: depends on prime)

    OUTPUT:

    a dictionary representing the mod `p` Adem relations
    applied to `P^a P^b` or (if `c` present) to `P^a \beta^b P^c`.

    EXAMPLES:

    If two arguments (`a` and `b`) are given, then computations are
    done mod 2.  If `a \geq 2b`, then the dictionary {(a,b): 1} is
    returned.  Otherwise, the right side of the mod 2 Adem relation
    for `\text{Sq}^a \text{Sq}^b` is returned.  For example, since
    `\text{Sq}^2 \text{Sq}^2 = \text{Sq}^3 \text{Sq}^1`, we have::

        sage: from sage.algebras.steenrod.steenrod_algebra_mult import adem
        sage: adem(2,2) # indirect doctest
        {(3, 1): 1}
        sage: adem(4,2)
        {(4, 2): 1}
        sage: adem(4,4) == {(6, 2): 1, (7, 1): 1}
        True

    If `p` is given and is odd, then with two inputs `a` and `b`, the
    Adem relation for `P^a P^b` is computed.  With three inputs `a`,
    `b`, `c`, the Adem relation for `P^a \beta^b P^c` is computed.
    In either case, the keys in the output are all tuples of odd length,
    with ``(i_1, i_2, ..., i_m)`` representing

    .. MATH::

        \beta^{i_1} P^{i_2} \beta^{i_3} P^{i_4} ... \beta^{i_m}

    For instance::

        sage: adem(3,1, p=3)
        {(0, 3, 0, 1, 0): 1}
        sage: adem(3,0,1, p=3)
        {(0, 3, 0, 1, 0): 1}
        sage: adem(1,0,1, p=7)
        {(0, 2, 0): 2}
        sage: adem(1,1,1, p=5) == {(0, 2, 1): 1, (1, 2, 0): 1}
        True
        sage: adem(1,1,2, p=5) == {(0, 3, 1): 1, (1, 3, 0): 2}
        True
    """
    if generic is None:
        generic = False if p==2 else True
    if not generic:
        return adem_2(a,b,c)
    else:
        return adem_odd(a, b, c, p)

@memoized
def make_mono_admissible_2(mono):
    if len(mono) == 1:
        return {mono: 1}
    if len(mono) == 2:
        return adem_2(*mono)
    # check to see if admissible:
    admissible = True
    for j in range(len(mono)-1):
        if mono[j] < 2*mono[j+1]:
            admissible = False
            break
    if admissible:
        return {mono: 1}
    # else j is the first index where admissibility fails
    ans = {}
    y = adem(mono[j], mono[j+1])
    for x in y:
        new = mono[:j] + x + mono[j+2:]
        new = make_mono_admissible_2(new)
        for m in new:
            if m in ans:
                ans[m] = ans[m] + y[x] * new[m]
                if ans[m] % p == 0:
                    del ans[m]
            else:
                ans[m] = y[x] * new[m]
    return ans

@memoized
def make_mono_admissible_odd(mono, p):
    # check to see if admissible:
    admissible = True
    for j in range(1, len(mono)-2, 2):
        if mono[j] < mono[j+1] + p*mono[j+2]:
            admissible = False
            break
    if admissible:
        return {mono: 1}
    # else j is the first index where admissibility fails
    ans = {}
    y = adem_odd(*mono[j:j+3], p=p)
    for x in y:
        new_x = list(x)
        new_x[0] = mono[j-1] + x[0]
        if len(mono) >= j+3:
            new_x[-1] = mono[j+3] + x[-1]
        if new_x[0] <= 1 and new_x[-1] <= 1:
            new = mono[:j-1] + tuple(new_x) + mono[j+4:]
            new = make_mono_admissible_odd(new, p)
            for m in new:
                if m in ans:
                    ans[m] += y[x] * new[m]
                    ans[m] = ans[m] % p
                    if ans[m] == 0:
                        del ans[m]
                else:
                    ans[m] = y[x] * new[m]
    return ans


def make_mono_admissible(mono, p=2, generic=None):
    r"""
    Given a tuple ``mono``, view it as a product of Steenrod
    operations, and return a dictionary giving data equivalent to
    writing that product as a linear combination of admissible
    monomials.

    When `p=2`, the sequence (and hence the corresponding monomial)
    `(i_1, i_2, ...)` is admissible if `i_j \geq 2 i_{j+1}` for all
    `j`.

    When `p` is odd, the sequence `(e_1, i_1, e_2, i_2, ...)` is
    admissible if `i_j \geq e_{j+1} + p i_{j+1}` for all `j`.

    INPUT:

    - ``mono`` - a tuple of non-negative integers
    - `p` - prime number, optional (default 2)
    - `generic` - whether to use the generic Steenrod algebra, (default: depends on prime)

    OUTPUT:

    Dictionary of terms of the form (tuple: coeff), where
    'tuple' is an admissible tuple of non-negative integers and
    'coeff' is its coefficient.  This corresponds to a linear
    combination of admissible monomials.  When `p` is odd, each tuple
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

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_mult import make_mono_admissible
        sage: make_mono_admissible((12,)) # already admissible, indirect doctest
        {(12,): 1}
        sage: make_mono_admissible((2,1)) # already admissible
        {(2, 1): 1}
        sage: make_mono_admissible((2,2))
        {(3, 1): 1}
        sage: make_mono_admissible((2, 2, 2))
        {(5, 1): 1}
        sage: make_mono_admissible((0, 2, 0, 1, 0), p=7)
        {(0, 3, 0): 3}

    Test the fix from :trac:`13796`::

        sage: SteenrodAlgebra(p=2, basis='adem').Q(2) * (Sq(6) * Sq(2)) # indirect doctest
        Sq^10 Sq^4 Sq^1 + Sq^10 Sq^5 + Sq^12 Sq^3 + Sq^13 Sq^2
    """
    if generic is None:
        generic = False if p==2 else True
    if len(mono) == 1:
        return {mono: 1}
    if not generic and len(mono) == 2:
        return adem(*mono, p=p, generic=generic)
    if not generic:
        return make_mono_admissible_2(mono)
    else:
        return make_mono_admissible_odd(mono, p)

def product_2(m1, m2):
    return make_mono_admissible_2(m1 + m2)

def product_odd(m1, m2, p):
    return make_mono_admissible_odd(m1 + m2, p)

def product(m1, m2, p, generic = None):
    return make_mono_admissible(m1 + m2, p, generic)




@memoized
def basis_2(n, bound = 1):
    # Build basis recursively.  last = last term.
    # last is >= bound, and we will append (last,) to the end of
    # elements from basis (n - last, bound=2 * last).
    # This means that 2 last <= n - last, or 3 last <= n.
    if(n == 0):
        return [[]]
    result = [[n]];
    for last in range(bound, 1 + n // 3):
        for vec in basis_2(n - last, 2 * last):
            result.append(vec + [last])
    return result
    
@memoized
def basis_odd(n, p, bound = 1):
    if(n == 0):
        return ((),)
    elif n % (2 * (p-1)) == 0 and n // (2 * (p-1)) >= bound:
        result = [(0, n // (2 * (p-1)), 0)]
    elif n == 1:
        result = [(1,)]
    else:
        result = []
    # 2 cases: append P^{last}, or append P^{last} beta
    # case 1: append P^{last}
    # We do decimal division and take ceiling so if (2*(p - 1)) divides n evenly, we get one fewer iteration.
    # Annoying that range cannot take floating endpoints. Note we need __future__ division for this to work.
    for last in range(bound, int(ceiling(n / (2*(p - 1))))):
        remaining_degree = n - 2*(p-1)*last
        basis = basis_odd(remaining_degree, p, p * last)
        for vec in basis:
            result.append(vec + (last, 0))
            
    # case 2: append P^{last} beta
    if(bound == 1):
        bound = 0
    # Note that in this loop we do integer division so we don't have the 
    # "it divides evenly" edge case that we had in the other one
    for last in range(bound + 1, 1 + (n // (2*(p - 1)))):
        remaining_degree = n - 2 * (p - 1) * last - 1
        basis = basis_odd(remaining_degree, p, p * last)
        for vec in basis:
            if not vec:
                vec = (0,)
            result.append(vec + (last, 1))
    return tuple(result)

def basis(n, p = 2, generic = None):
    """
    Serre-Cartan basis in dimension `n`.

    INPUT:

    - ``n`` - non-negative integer
    - ``bound`` - positive integer (optional)
    - ``prime`` - positive prime number (optional, default 2)

    OUTPUT: tuple of mod p Serre-Cartan basis elements in dimension n

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_bases import serre_cartan_basis
        sage: serre_cartan_basis(7)
        ((7,), (6, 1), (4, 2, 1), (5, 2))
        sage: serre_cartan_basis(13,3)
        ((1, 3, 0), (0, 3, 1))
        sage: serre_cartan_basis(50,5)
        ((1, 5, 0, 1, 1), (1, 6, 1))
    """
    generic = generic or p != 2;
    if not generic: 
        return basis_2(n)
    else:
        return basis_odd(n, p)
