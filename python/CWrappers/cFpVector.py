from ctypes import *
from ctypes_wrap import *

class cVector:
    def __init__(self, p=None, dim=None, offset=0, vector=None):
        self.freed = False
        if vector != None and type(vector) != list:
            self.p = vector.contents.p
            self.dimension = vector.contents.dimension
            self.c_list_type = c_uint * self.dimension
            self.vector = vector
            return
        
        if type(vector) == list:
            dim = len(vector)

        self.p = p
        self.dimension = dim
        self.c_list_type = c_uint * self.dimension
        CSteenrod.initializePrime(p)
        self.vector = CSteenrod.Vector_construct(p, dim, offset)

        if type(vector) == list:
            self.pack(vector)      

    def free(self):
        if not self.freed:
            CSteenrod.Vector_free(self.vector)
        self.freed = True

    def assign(self, w):
        CSteenrod.Vector_assign(self.vector, w.vector)
    
    def pack(self, list):
        if self.dimension != len(list):
            raise Exception("Wrong length.")
        c_list = self.c_list_type()
        for i, elt in enumerate(list):
            c_list[i] = elt % self.p
        CSteenrod.Vector_pack(self.vector, c_list)

    def unpack(self):
        c_list = self.c_list_type()
        CSteenrod.Vector_unpack(c_list, self.vector)
        py_list = [None] * self.dimension
        for i in range(self.dimension):
            py_list[i] = c_list[i]
        return py_list        
    
    def addBasisElement(self, idx, c=1):
        c = c % self.p
        CSteenrod.Vector_addBasisElement(self.vector, idx, c)

    def add(self, w, c=1):
        c = c % self.p
        CSteenrod.Vector_add(self.vector, w.vector, c)
        
        
    def scale(self, c):
        c = c % self.p
        CSteenrod.Vector_scale(self.vector, c)

    def slice(self, min, max):
        cSlice = CSteenrod.Vector_construct(self.p, 0, 0)
        CSteenrod.Vector_slice(cSlice, self.vector, min, max)
        slice = cVector(vector=cSlice)
        return slice

    def __len__(self):
        return self.dimension
        
    def __setitem__(self, idx, value):
        CSteenrod.Vector_setEntry(self.vector, idx, (value % self.p))

    def __getitem__(self, key):
        if isinstance(key, slice):
            #Get the start, stop, and step from the slice
            if(key.step != None and key.step != 1):
                print(key.step)
                raise(IndexError("Slice steps not equal to 1 not supported."))
            return self.slice(key.start, key.stop)
        elif isinstance(key, int):
            if key < 0 : #Handle negative indices
                key += len(self)
            if key < 0 or key >= len(self):
                raise(IndexError("The index (%d) is out of range."%key))
            return CSteenrod.Vector_getEntry(self.vector, key)
        else:
            raise TypeError("Invalid argument type.")

    def __iter__(self):
        return cVector_iterator(self)


class cVector_iterator:
    def __init__(self, vector):
        self.cIterator = CSteenrod.Vector_getIterator(vector.vector)

    def __next__(self):
        if not self.cIterator.has_more:
            raise StopIteration
        result = self.cIterator.value
        cIterator = CSteenrod.Vector_stepIterator(self.cIterator)
        return result


class cMatrix:
    def __init__(self, p, rows, columns):
        self.p = p
        self.rows = rows
        self.columns = columns
        CSteenrod.initializePrime(p)
        self.cM = CSteenrod.Matrix_construct(p, rows, columns)
    
    def pack(self, py_M):
        for i in range(c_M.contents.rows):
            c_packVector(c_M.contents.matrix[i], py_M[i])

    def unpack(c_M):
        return [c_unpackVector(c_M.contents.matrix[i]) for i in range(c_M.contents.rows)]


def c_row_reduce(c_M):
    array_type = c_int * c_M.contents.columns
    pivots_array = array_type()
    CSteenrod.rowReduce(c_M, pivots_array, c_M.contents.rows)
    c_M.pivots = pivots_array

def vector_to_C(p, vector):
    c_v = cVector(p, len(vector))
    cVector_pack(c_v, vector)
    return c_v

def matrix_to_C(p, matrix):
    rows = len(matrix)
    columns = len(matrix[0])
    c_M = cMatrix_construct(p, rows, columns)
    cMatrix_pack(c_M, matrix)
    return c_M

def matrix_from_C(matrix):
    return cMatrix_unpack(matrix)


if __name__ == "__main__":
    p = 2
    dim = 14

    
    
    

