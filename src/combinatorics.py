def base_p_expansion(n, p, padlength = 0):
    result = []
    while n != 0:
        result.append(n % p)
        n /= p
    while len(result) < padlength:
        result.append(0)
    return result

# Integral binomial coefficient.
def binomial(n, k):
    coeff = reduce( lambda x, y: x*y, range(n-k+1, n+1), 1 )
    for x in range(1,k+1):
        coeff /= x
    return coeff


# Mod 2 multinomial coefficient
def multinomial_mod_2(list):
    old_sum = list[0]
    for x in list[1:]:
        j = 1
        while j <= min(old_sum, list[i]):
            if ((j & old_sum ) == j) and ((j & list[i]) != 0):
                return 0
            j = j << 1
        old_sum += list[i]
    return 1

def binomial_mod2(n,k):
    if n < k:
        return 0
    else:
        return + ( (n-k) & k == 0 )


def multinomial_odd(list,p):
    n = sum(list)
    answer = 1
    n_expansion = base_p_expansion(n, p)
    list_expansion = map(lambda x: base_p_expansion(x, p, len(n_expansion) ), list)
    index = 0
    for index in range(len(n_expansion)):
        multi = 1
        partial_sum = 0
        for exp in list_expansion:
            if index < len(exp):
                partial_sum += (+exp[index])
                multi *= binomial(partial_sum, +exp[index])
        answer = (answer * multi) % p
        if answer == 0:
            return 0
    return answer

def binomial_modp(n, k, p):
    if n < k:
        return 0
    return multinomial_odd([n-k, k], p)


def multinomial(list, p = 2):
    if p == 2:
        return multinomial_mod_2(list)
    else:
        return multinomial_odd(list, p)


def binomial_gen(n, k, p = 2):
    if n<k or k<0:
        return 0
    if p == 2:
        return binomial_mod2(n,k)
    else:
        return binomial_modp(n,k,p)



def restricted_partitions(n, l, no_repeats = False):
    """
            restricted_partitions(10, [6,4,2])
            [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
            restricted_partitions(10, [6,4,2,2,2])
            [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
            restricted_partitions(10, [6,4,4,4,2,2,2,2,2,2])
            [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
    """
    if n < 0:
        return []
    elif n == 0:
        return [[]]
    else:
        results = []
        index = +no_repeats
        old_i = 0
        for i in l:
            if old_i != i:
                for sigma in restricted_partitions(n-i, l[ index : ], no_repeats):
                    results.append([i] +  sigma)
            index += 1
            old_i = i
        return results

def xi_degrees(n, p=2, reverse = True):
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
    xi_max = 1
    while N > 0:
        N = N/p
        xi_max += 1
    l = []
    for d in range(1, xi_max):
        l.append(((p**d-1)/(p-1)))
    if(reverse):
        l.reverse()
    return l


def WeightedIntegerVectors(n, l):
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
        
    if len(l) == 1:
        if n % l[0] == 0:
            yield [n / l[0]]
        return

    k = 0
    cur = [ n / l[k]  + 1 ]
    rem = n - cur[ -1 ] * l[k] # Amount remaining
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
                yield cur + [rem / l[- 1]]
        else:
            k += 1
            cur.append(rem / l[k] + 1)
            rem -= cur[- 1] * l[k]

# fle.log(Array.from(WeightedIntegerVectors(3, [2,1,1])))
