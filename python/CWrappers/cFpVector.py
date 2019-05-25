from ctypes import *
from ctypes_wrap import *

def construct_c_vector(p, dim):
    CSteenrod.initializePrime(p)
    v = CSteenrod.constructVector(p, dim)
    v.p = p
    v.dimension = dim
    v.freed = False
    return v
    
def free_c_vector(v):
    if type(v) == c_vector and not v.freed:
        CSteenrod.freeVector(v)
        
def assign_c_vector(v, w):
    CSteenrod.assignVector(v.p, v, w)
    
def c_packVector(v, list):
    if v.dimension != len(list):
        raise Exception("Wrong length.")
    c_list_type = c_uint * len(list)
    c_list = c_list_type()
    for i, elt in enumerate(list):
        c_list[i] = elt % v.p
    CSteenrod.packVector(v.p, v, c_list)

def c_unpackVector(v):
    c_list_type = c_uint * v.dimension
    c_list = c_list_type()
    CSteenrod.unpackVector(v.p, c_list, v)
    py_list = [None] * v.dimension
    for i in range(v.dimension):
        py_list[i] = c_list[i]
    return py_list
    
def c_getVectorEntry(v, idx):
    return CSteenrod.getVectorEntry(v.p, v, idx)

def c_setVectorEntry(v, idx, value):
    CSteenrod.setVectorEntry(v.p, v, idx, (value % v.p))
    
def c_addBasisElementToVector(v, idx, c=1):
    c = c % v.p
    if v.p == 2:
        CSteenrod.addBasisElementToVector2(v.p, v, idx, c)
    else:
        CSteenrod.addBasisElementToVectorGeneric(v.p, v, idx, c)    

def c_addVectors(v, w, c=1):
    c = c % v.p
    if v.p == 2:
        CSteenrod.addVectors2(v.p, v, w, c)
    else:
        CSteenrod.addVectorsGeneric(v.p, v, w, c)        
    
    
def c_scaleVector(v, c):
    c = c % v.p
    if v.p == 2:
        CSteenrod.scaleVector2(v.p, v, c)
    else:
        CSteenrod.scaleVectorGeneric(v.p, v, c)        
    
    
if __name__ == "__main__":
    import random
    p = 2
    dim = 14
    v = construct_c_vector(5, dim)
    w = construct_c_vector(5, dim)
    k = [random.randint(0,p-1) for x in range(dim)]
    l = [random.randint(0,p-1) for x in range(dim)]
    result = [ (k[i] + l[i]) % 5 for i in range(dim)]
    c_packVector(v, k)
    c_packVector(w, l)
    c_addVectors(v, w)
    s = c_unpackVector(v)
    for x in range(dim):
        if c_getVectorEntry(v, x) != s[x]:
            print("getVectorEntry maybe failed")
            break
        
    if result != s:
        print("Test failed")
    
    
    

