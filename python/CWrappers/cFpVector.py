from ctypes import *
from ctypes_wrap import *

def construct_c_vector(p, dim, offset=0):
    CSteenrod.initializePrime(p)
    if(p==2):
        v = CSteenrod.constructVector2(p, dim, offset)
    else:
        v = CSteenrod.constructVectorGeneric(p, dim, offset)
    v.p = p
    v.dimension = dim
    v.freed = False
    return v
    
def free_c_vector(v):
    if type(v) == c_vector and not v.freed:
        CSteenrod.freeVector(v)
    v.freed = True
        
def assign_c_vector(v, w):
    CSteenrod.assignVector(v, w)
    
def c_packVector(v, list):
    if v.dimension != len(list):
        raise Exception("Wrong length.")
    c_list_type = c_uint * len(list)
    c_list = c_list_type()
    for i, elt in enumerate(list):
        c_list[i] = elt % v.p
    CSteenrod.packVector(v, c_list)

def c_unpackVector(v):
    c_list_type = c_uint * v.dimension
    c_list = c_list_type()
    CSteenrod.unpackVector(c_list, v)
    py_list = [None] * v.dimension
    for i in range(v.dimension):
        py_list[i] = c_list[i]
    return py_list
    
def c_getVectorEntry(v, idx):
    return CSteenrod.getVectorEntry(v, idx)

def c_setVectorEntry(v, idx, value):
    CSteenrod.setVectorEntry(v, idx, (value % v.p))
    
def c_addBasisElementToVector(v, idx, c=1):
    c = c % v.p
    if v.p == 2:
        CSteenrod.addBasisElementToVector2(v, idx, c)
    else:
        CSteenrod.addBasisElementToVectorGeneric(v, idx, c)    

def c_addVectors(v, w, c=1):
    c = c % v.p
    if v.p == 2:
        CSteenrod.addVectors2(v, w, c)
    else:
        CSteenrod.addVectorsGeneric(v, w, c)        
    
    
def c_scaleVector(v, c):
    c = c % v.p
    if v.p == 2:
        CSteenrod.scaleVector2(v, c)
    else:
        CSteenrod.scaleVectorGeneric(v, c)        
    
    
def test_c_vector(p, dim):
    import random
    v = construct_c_vector(p, dim)
    w = construct_c_vector(p, dim)
    k = [random.randint(0,p-1) for x in range(dim)]
    l = [random.randint(0,p-1) for x in range(dim)]
    print(k)
    print(l)
    result = [ (k[i] + l[i]) % p for i in range(dim)]
    print(result)
    c_packVector(v, k)
    k_packed_unpacked = c_unpackVector(v)
    if k != k_packed_unpacked:
        print("Pack unpack failed.")
        print("Orig: ", k)
        print("Punp: ", k_packed_unpacked)
        return

    c_packVector(w, l)
    c_addVectors(v, w)
    s = c_unpackVector(v)
    for x in range(dim):
        if c_getVectorEntry(v, x) != s[x]:
            print("getVectorEntry and unpack disagree")
            break
    if result != s:
        print("Test failed:")
        print("k:      ", k)
        print("l:      ", l)
        print("result: ", result)
        print("s:      ", s)


if __name__ == "__main__":
    p = 2
    dim = 14

    
    
    

