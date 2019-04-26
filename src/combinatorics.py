from memoized import memoized

def base_p_expansion(n, p, padlength = 0):
    """Return the base p expansion of n as a list [d_0, d_1, d_2, ...] so that sum(d_i p^i) = n"""
    result = []
    while n != 0:
        result.append(n % p)
        n //= p
    while len(result) < padlength:
        result.append(0)
    return result

binomial_table = {}

def direct_binomial_initialize_table(p):
    """Produce the table with n choose k mod p for n, k < p"""
    table_p = [[0] * p for _ in range(p)]
    for n in range(p):
        entry = 1
        table_p[n][0] = entry
        for k in range(1, n + 1):
            entry *= (n + 1 - k)
            entry //= k
            table_p[n][k] = entry % p
    binomial_table[p] = table_p

def direct_binomial(n, k, p):
    """Direct binomial coefficient, used for n, k < p"""
    if p not in binomial_table:
        direct_binomial_initialize_table(p)
    return binomial_table[p][n][k]

# Mod 2 multinomial coefficient
#@memoized # No point in memoizing unless we're passing it tuples
def multinomial_2(l):
    """The mod 2 multinomial coefficient. Implemented as Sum(l) = BitOr(l)."""
    bit_or = 0
    sum = 0
    for v in l:
        sum += v
        bit_or |= v
        if bit_or < sum:
            return 0           
    return 1

@memoized
def binomial_2(n,k):
    """The mod 2 binomial coefficient. If n - k and k have no common 1's in binary, return 1 else return 0"""
    if n < k:
        return 0
    else:
        return  +((n-k) & k == 0)


def multinomial_odd(l, p):
    """The multinomial coefficient of l at an odd prime."""
    n = sum(l)
    answer = 1
    n_expansion = base_p_expansion(n, p)
    l_expansions = [base_p_expansion(x, p, len(n_expansion)) for x in l]
    index = 0
    for index in range(len(n_expansion)):
        multi = 1
        partial_sum = 0
        for exp in l_expansions:
            partial_sum += exp[index]
            if partial_sum > n_expansion[index]:
                return 0
            multi *= direct_binomial(partial_sum, exp[index], p)
        answer = (answer * multi) % p
    return answer

def binomial_odd(n, k, p):
    """Binomial odd: just use multinomial_odd on [n-k, k]."""
    if n < k:
        return 0
    return multinomial_odd([n-k, k], p)


def multinomial(l, p): 
    """Dispatch to multinomial_2 or multinomial_odd."""
    if p == 2:
        return multinomial_2(l)
    else:
        return multinomial_odd(l, p)


def binomial(n, k, p):
    """Dispatch to binomial_2 or binomial_odd."""  
    if n<k or k<0:
        return 0
    if p == 2:
        return binomial_2(n,k)
    else:
        return binomial_odd(n,k,p)


def xi_degrees(n, *, p):
    """Return the degrees of the xi_i's in dimension <= n in increasing order.
       xi_degrees(17, p=2) --> [1, 3, 7, 15]
       xi_degrees(17, p=3) --> [1, 4, 13]
       xi_degrees(400, p=17) --> [1, 18, 307]
    """
    if n <= 0:
        return []
    N = n*(p-1) + 1
    xi_max = 0
    while N > 0:
        N //= p
        xi_max += 1
    result = [0] * (xi_max - 1)
    entry = 0
    p_to_the_d = 1
    for i in range(xi_max - 1):
        entry += p_to_the_d
        p_to_the_d *= p    
        result[i] = entry
    return result
    
def tau_degrees(n, *, p):
    """Return the degrees of tau_i's in dimension <= n in increasing order.
       If p == 2, this has the same output as xi_degrees.
    """
    xi_degs = [0] + xi_degrees((n - 1)//(2*(p-1)), p=p)
    return [1+2*(p-1)*d for d in xi_degs]

def min_if_b_not_zero_else_a(a, b):
    """This is a helper function for WeightedIntegerVectors. 
       Since weights must be positive, we use 0 for max_weight to denote no max weight.
       Thus, we need a min function that's a no-op when the second argument is 0.
       It also returns True if the result is equal to the first argument, False if not
    """
    
    if a <= b or b == 0:
        return (a, True)
    return (b, False)

def WeightedIntegerVectors(n, l, max_weight = 0):
    """
    Iterate over all ``l`` weighted integer vectors with total weight ``n``.
    Warning: it does this in place.

    INPUT:

    - ``n`` -- an integer
    - ``l`` -- the weights in weakly increasing order

    EXAMPLES::

          list(WeightedIntegerVectors(3, [1, 1,  2]))
          --> [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]] 
        
        Oops! Let's try again....
        
          [tuple(x) for x in WeightedIntegerVectors(3, [1, 1, 2])]
          --> [(0, 1, 1), (1, 0, 1), (0, 3, 0), (1, 2, 0), (2, 1, 0), (3, 0, 0)]
        
        That's more like it.
    """
    if n < 0:
        return

    if not l:
        if n == 0:
            yield []
        return
        
    # We might be in case 3 below right away. Cover it separately up here.
    if len(l) == 1:
        rem = n
        if rem % l[0] == 0:
            ratio, ratio_leq_max_weight = min_if_b_not_zero_else_a(rem // l[0], max_weight)
            if ratio_leq_max_weight:
                yield [ratio]
        return

    cur = [0] * len(l)
    k = len(l) - 1
    # Now we're definitely in "case 4"
    rem = n
    ratio, _ = min_if_b_not_zero_else_a( rem // l[k], max_weight)
    cur[k] = ratio + 1
    rem -= cur[ k ] * l[k] # Amount remaining
    while k < len(cur):
        cur[k] -= 1
        rem += l[k]
        if rem == 0:
            yield cur
        elif cur[k] < 0 or rem < 0 :
            rem += cur[k] * l[k]
            cur[k] = 0
            k += 1
        elif k == 1:
            if rem % l[0] == 0:
                ratio, ratio_leq_max_weight = min_if_b_not_zero_else_a(rem // l[0], max_weight)
                if ratio_leq_max_weight:
                    cur[0] = ratio
                    yield cur
                    cur[0] = 0
        else:
            k -= 1
            ratio, _ = min_if_b_not_zero_else_a(rem // l[k], max_weight)
            cur[k] = ratio + 1
            rem -= cur[k] * l[k]


def restricted_partitions(n, l):
    """
            restricted_partitions(10, [6,4,2])
            [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
            restricted_partitions(10, [6,4,2,2,2])
            [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
            restricted_partitions(10, [6,4,4,4,2,2,2,2,2,2])
            [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
    """
    return WeightedIntegerVectors(n, l, 1)
        
