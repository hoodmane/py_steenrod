import pytest
import random

from cFpVector import *

primes = [2,3,5,7]
dimensions = [5, 10, 33, 65, 1000]
repeats = 100

@pytest.mark.parametrize("dim", dimensions)
@pytest.mark.parametrize("p", primes)
def test_pack_unpack(p, dim):
    v = cVector(p, dim)
    k = [random.randint(0,p-1) for x in range(dim)]
    v.pack(k)
    k_packed_unpacked = v.unpack()
    assert k == k_packed_unpacked
    v.free()

@pytest.mark.parametrize("dim", dimensions)
@pytest.mark.parametrize("p", primes)
def test_pack_get(p, dim):
    v = cVector(p, dim)
    k = [random.randint(0,p-1) for x in range(dim)]
    v.pack(k)
    k_packed_unpacked = v.unpack()
    for i in range(dim):
        assert v[i] == k[i]
    v.free()


@pytest.mark.parametrize("dim", dimensions)
@pytest.mark.parametrize("p", primes)
def test_set_get(p, dim):
    v = cVector(p, dim)
    k = [0] * dim
    for i in range(repeats):
        index = random.randint(0, dim-1)
        value = random.randint(0, p-1)
        k[index] = value
        v[index] = value
        for j in range(dim):
            assert v[index] == k[index]
    result = v.unpack()
    v.free()
    assert result == k


@pytest.mark.parametrize("dim", dimensions)
@pytest.mark.parametrize("p", primes)
def test_assign(p, dim):
    v = cVector(p, dim)
    w = cVector(p, dim)
    k = [random.randint(0,p-1) for x in range(dim)]
    l = [random.randint(0,p-1) for x in range(dim)]
    v.pack(k)
    w.pack(l)
    v.assign(w)
    result = v.unpack()
    v.free()
    w.free()
    assert result == l

@pytest.mark.parametrize("dim", dimensions)
@pytest.mark.parametrize("p", primes)
def test_addBasisElement(p, dim):
    v = cVector(p, dim)
    k = [0] * dim
    for i in range(repeats):
        index = random.randint(0, dim-1)
        value = random.randint(0, p-1)
        k[index] += value
        k[index] = k[index] % p
        v.addBasisElement(index, value)
    result = v.unpack()
    v.free()    
    assert result == k


@pytest.mark.parametrize("dim", dimensions)
@pytest.mark.parametrize("p", primes)
def test_add(p, dim):
    v = cVector(p, dim)
    w = cVector(p, dim)
    k = [random.randint(0,p-1) for x in range(dim)]
    l = [random.randint(0,p-1) for x in range(dim)]
    py_sum = [ (k[i] + l[i]) % p for i in range(dim)]
    cVector.pack(v, k)
    cVector.pack(w, l)
    cVector.add(v, w)
    result = cVector.unpack(v)
    assert result == py_sum
    v.free()
    w.free()

@pytest.mark.parametrize("dim", dimensions)
@pytest.mark.parametrize("p", primes)
def test_addArray(p, dim):
    v = cVector(p, dim)
    w = (c_uint * dim)()
    k = [random.randint(0,p-1) for i in range(dim)]
    l = [random.randint(0,p-1) for i in range(dim)]
    py_sum = [ (k[i] + l[i]) % p for i in range(dim)]
    for i in range(dim):
        w[i] = l[i]
    cVector.pack(v, k)
    CSteenrod.Vector_addArray(v.vector, w, 1)
    result = cVector.unpack(v)
    assert result == py_sum
    v.free()

@pytest.mark.parametrize("dim,min,max", [(10, 3,7), (100, 57, 72), (1000, 178, 256), (1000, 700, 900)])
@pytest.mark.parametrize("p", primes)
def test_slice(p, dim, min, max):
    v = cVector(p, dim)
    k = [random.randint(0,p-1) for x in range(dim)]
    v.pack(k)
    print("slicing")
    w = v[min:max]
    print("sliced")
    result = w.unpack()
    print("seg?")
    assert result == k[min:max]
    v.free()
    w.free()

