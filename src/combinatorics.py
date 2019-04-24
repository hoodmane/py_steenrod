from memoized import memoized

def base_p_expansion(n, p, padlength = 0):
    result = []
    while n != 0:
        result.append(n % p)
        n //= p
    while len(result) < padlength:
        result.append(0)
    return result

binomial_table = {}

def direct_binomial_initialize_table(p):
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
    """direct binomial coefficient, used for n, k < p"""
    if p not in binomial_table:
        direct_binomial_initialize_table(p)
    table = binomial_table[p]
    return table[n][k]

# Mod 2 multinomial coefficient
#@memoized # No point in memoizing unless we're passing it tuples
def multinomial_2(l):
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
    if n < k:
        return 0
    else:
        return  +((n-k) & k == 0)


def multinomial_odd(l, p):
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
    if n < k:
        return 0
    return multinomial_odd([n-k, k], p)


def multinomial(l, p):
    if p == 2:
        return multinomial_2(l)
    else:
        return multinomial_odd(l, p)


def binomial(n, k, p):
    if n<k or k<0:
        return 0
    if p == 2:
        return binomial_2(n,k)
    else:
        return binomial_odd(n,k,p)


def xi_degrees(n,*, p, reverse):
    """
            sage: sage.algebras.steenrod.steenrod_algebra_bases.xi_degrees(17)
            [15, 7, 3, 1]
            sage: sage.algebras.steenrod.steenrod_algebra_bases.xi_degrees(17, reverse=False)
            [1, 3, 7, 15]
            sage: sage.algebras.steenrod.steenrod_algebra_bases.xi_degrees(17,p=3)
            [13, 4, 1]
            sage: sage.algebras.steenrod.steenrod_algebra_bases.xi_degrees(400,p=17)
            [307, 18, 1]
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
    if(reverse):
        result.reverse()
    return result


def min_if_b_not_zero_else_a(a, b):
    if a <= b or b == 0:
        return (a, True)
    return (b, False)

def WeightedIntegerVectors(n, l, max_weight = 0):
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
    if n < 0:
        return

    if not l:
        if n == 0:
            yield []
        return
        
    # We might be in case 3 below right away. Cover it separately up here.
    if len(l) == 1:
        rem = n
        if rem % l[-1] == 0:
            _, ratio_leq_max_weight = min_if_b_not_zero_else_a(rem // l[-1], max_weight)
            if ratio_leq_max_weight:
                yield [rem // l[-1]]
        return

    k = 0
    # Now we're definitely in "case 4"
    rem = n
    ratio, _ = min_if_b_not_zero_else_a( rem // l[k], max_weight)
    cur = [ ratio + 1 ]
    rem -= cur[ -1 ] * l[k] # Amount remaining
    while len(cur) > 0:
        cur[-1] -= 1
        rem += l[k]
        if rem == 0:
            yield cur + [0 for _ in range(len(l) - len(cur))]
        elif cur[- 1] < 0 or rem < 0 :
            rem += cur.pop() * l[k]
            k -= 1
        elif len(l) == len(cur) + 1:
            if rem % l[-1] == 0:
                ratio, ratio_leq_max_weight = min_if_b_not_zero_else_a(rem // l[-1], max_weight)
                if ratio_leq_max_weight:
                    yield cur + [ratio]
        else:
            k += 1
            ratio, _ = min_if_b_not_zero_else_a(rem // l[k], max_weight)
            cur.append(ratio + 1)
            rem -= cur[-1] * l[k]


def restricted_partitions(n, l):
    """
            restricted_partitions(10, [6,4,2])
            [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
            restricted_partitions(10, [6,4,2,2,2])
            [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
            restricted_partitions(10, [6,4,4,4,2,2,2,2,2,2])
            [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
    """
    for weights in WeightedIntegerVectors(n, l, 1):
        result = []
        for idx, i in enumerate(weights):
            if i == 1:
                result.append(l[idx])
        yield result
        
