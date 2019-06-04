from ctypes import *

# typedef struct  {
#     uint dimension;
#     uint size;
#     uint offset;
# } Vector;
class c_Vector(Structure):
    _fields_ = [
        ("dimension", c_uint),
        ("size", c_uint),
        ("offset", c_uint)
    ]

# typedef struct {
#     bool has_more;
#     uint index;
#     uint value;
# // private
#     Vector * vector;
#     uint limb_index;
#     uint bit_index;    
# } VectorIterator;

class c_VectorIterator(Structure):
    _fields_ = [
        ("has_more",c_bool),
        ("index", c_uint),
        ("value", c_uint),
        ("_vector", c_Vector),
        ("_limb_index", c_uint),
        ("_bit_index", c_uint)
    ]

# typedef struct {
#     uint p;
#     uint rows;
#     uint columns;
#     Vector ** matrix;
# } Matrix;
class c_Matrix(Structure):
    _fields_ = [
        ("p", c_uint),
        ("rows", c_uint),
        ("columns", c_uint),
        ("matrix", POINTER(POINTER(c_Vector)))
    ]

def wrap_FpVector(CSteenrod):        
    # Vector * Vector_initialize(uint p, uint64 * vector_container, uint64 * memory, uint dimension, uint offset);
    # Vector * Vector_construct(uint p, uint dimension, uint offset);
    # void Vector_free(Vector * v); 

    CSteenrod.Vector_construct.argtypes = [c_uint, c_uint, c_uint]    
    CSteenrod.Vector_construct.restype = POINTER(c_Vector)

    CSteenrod.Vector_free.argtypes = [POINTER(c_Vector)]
    
    #void Vector_assign(Vector * target, Vector * source);
    CSteenrod.Vector_assign.argtypes = [POINTER(c_Vector), POINTER(c_Vector)]
    
    #void Vector_setToZero(Vector * target);
    CSteenrod.Vector_setToZero.argtypes = [POINTER(c_Vector)]

    #void Vector_pack(Vector * target, uint * source);
    CSteenrod.Vector_pack.argtypes = [POINTER(c_Vector), POINTER(c_uint)]
    
    #void Vector_unpack(uint * target, Vector * source);
    CSteenrod.Vector_unpack.argtypes = [POINTER(c_uint), POINTER(c_Vector)]
    
    #uint Vector_getEntry(Vector * v, uint index);
    CSteenrod.Vector_getEntry.argtypes = [POINTER(c_Vector), c_uint]    
    CSteenrod.Vector_getEntry.restype = c_uint
    
    #void Vector_setEntry(Vector * v, uint index, uint value);
    CSteenrod.Vector_setEntry.argtypes = [POINTER(c_Vector), c_uint, c_uint]    
    
    # void Vector_slice(Vector *result, Vector *source, uint start, uint end);
    # VectorIterator Vector_getIterator(Vector * v); 
    # VectorIterator Vector_stepIterator(VectorIterator);
    CSteenrod.Vector_slice.argtypes = [POINTER(c_Vector), POINTER(c_Vector), c_uint, c_uint]

    CSteenrod.Vector_getIterator.argtypes = [POINTER(c_Vector)]
    CSteenrod.Vector_getIterator.restype = c_VectorIterator
    CSteenrod.Vector_stepIterator.argtypes = [c_VectorIterator]
    CSteenrod.Vector_stepIterator.restype = c_VectorIterator

    #void Vector_addBasisElement(Vector * v, uint idx, uint c);
    CSteenrod.Vector_addBasisElement.argtypes = [POINTER(c_Vector), c_uint, c_uint]
    
    #void Vector_add(Vector * target, Vector * source, uint c);
    CSteenrod.Vector_add.argtypes = [POINTER(c_Vector), POINTER(c_Vector), c_uint]  
    
    #void Vector_addArray(Vector * target, uint *source, uint c);
    CSteenrod.Vector_addArray.argtypes = [POINTER(c_Vector), POINTER(c_uint), c_uint]

    #void Vector_scale(Vector * v, uint c);
    CSteenrod.Vector_scale.argtypes = [POINTER(c_Vector), c_uint] 
    
    
    #uint vectorToString(char * buffer, Vector * v);
    CSteenrod.Vector_toString.argtypes = [c_char_p, POINTER(c_Vector)]
    CSteenrod.Vector_toString.restype = c_uint
    
    #Vector* Matrix_construct(uint p, uint rows, uint cols);
    CSteenrod.Matrix_construct.argtypes = [c_uint, c_uint, c_uint]
    CSteenrod.Matrix_construct.restype = POINTER(c_Matrix)

    #void Matrix_free(Matrix *matrix);
    CSteenrod.Matrix_free.argtypes = [POINTER(c_Matrix)]

    #void rowReduce(Vector **matrix, int * column_to_pivot_row, uint rows);
    CSteenrod.rowReduce.argtypes = [POINTER(c_Matrix), POINTER(c_int)]