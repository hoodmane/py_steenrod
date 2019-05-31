from ctypes import *
from ctypes_wrap import *

def cVector_construct(p, dim, offset=0):
    CSteenrod.initializePrime(p)
    if(p==2):
        v = CSteenrod.Vector2_construct(p, dim, offset)
    else:
        v = CSteenrod.VectorGeneric_construct(p, dim, offset)
    v.freed = False
    return v
    
def cVector_free(v):
    if type(v) == c_Vector and not v.freed:
        CSteenrod.Vector_free(v)
    v.freed = True
        
def cVector_assign(v, w):
    CSteenrod.Vector_assign(v, w)
    
def cVector_pack(v, list):
    if v.contents.dimension != len(list):
        raise Exception("Wrong length.")
    c_list_type = c_uint * len(list)
    c_list = c_list_type()
    for i, elt in enumerate(list):
        c_list[i] = elt % v.contents.p
    CSteenrod.Vector_pack(v, c_list)

def cVector_unpack(v):
    c_list_type = c_uint * v.contents.dimension
    c_list = c_list_type()
    CSteenrod.Vector_unpack(c_list, v)
    py_list = [None] * v.contents.dimension
    for i in range(v.contents.dimension):
        py_list[i] = c_list[i]
    return py_list
    
def cVector_getEntry(v, idx):
    return CSteenrod.Vector_getEntry(v, idx)

def cVector_setEntry(v, idx, value):
    CSteenrod.Vector_setEntry(v, idx, (value % v.contents.p))
    
def cVector_addBasisElement(v, idx, c=1):
    c = c % v.p
    if v.contents.p == 2:
        CSteenrod.Vector2_addBasisElement(v, idx, c)
    else:
        CSteenrod.VectorGeneric_addBasisElement(v, idx, c)    

def cVector_add(v, w, c=1):
    c = c % v.contents.p
    if v.contents.p == 2:
        CSteenrod.Vector2_add(v, w, c)
    else:
        CSteenrod.VectorGeneric_addVectors(v, w, c)        
    
    
def cVector_scale(v, c):
    c = c % v.contents.p
    if v.contents.p == 2:
        CSteenrod.Vector2_scale(v, c)
    else:
        CSteenrod.VectorGeneric_scaleVector(v, c)        

def cMatrix_construct(p, rows, columns):
    CSteenrod.initializePrime(p)
    if p == 2:
        M = CSteenrod.Matrix2_construct(p, rows, columns)
    else:
        M = CSteenrod.MatrixGeneric_construct(p, rows, columns)
    return M

def cMatrix_pack(c_M, py_M):
    for i in range(c_M.contents.rows):
        c_packVector(c_M.contents.matrix[i], py_M[i])

def cMatrix_unpack(c_M):
    return [c_unpackVector(c_M.contents.matrix[i]) for i in range(c_M.contents.rows)]

def c_row_reduce(c_M):
    array_type = c_int * c_M.contents.columns
    pivots_array = array_type()
    CSteenrod.rowReduce(c_M, pivots_array, c_M.contents.rows)
    c_M.pivots = pivots_array

def vector_to_C(p, vector):
    c_v = cVector_construct(p, len(vector))
    cVector_pack(c_v, vector)
    return c_v

def vector_from_C(vector):
    return cVector_unpack(vector)

def matrix_to_C(p, matrix):
    rows = len(matrix)
    columns = len(matrix[0])
    c_M = cMatrix_construct(p, rows, columns)
    cMatrix_pack(c_M, matrix)
    return c_M

def matrix_from_C(matrix):
    return cMatrix_unpack(matrix)

    
def test_c_vector(p, dim):
    import random
    v = cVector_construct(p, dim)
    w = cVector_construct(p, dim)
    k = [random.randint(0,p-1) for x in range(dim)]
    l = [random.randint(0,p-1) for x in range(dim)]
    print(k)
    print(l)
    result = [ (k[i] + l[i]) % p for i in range(dim)]
    print(result)
    cVector_pack(v, k)
    k_packed_unpacked = cVector_unpack(v)
    if k != k_packed_unpacked:
        print("Pack unpack failed.")
        print("Orig: ", k)
        print("Punp: ", k_packed_unpacked)
        return

    cVector_pack(w, l)
    cVector_add(v, w)
    s = cVector_unpack(v)
    for x in range(dim):
        if cVector_getEntry(v, x) != s[x]:
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

    
    
    

