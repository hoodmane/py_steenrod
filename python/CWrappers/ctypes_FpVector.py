from ctypes import *

# typedef struct  {
#     struct VectorInterface_s * interface;
#     uint p;
#     uint dimension;
#     uint size;
#     uint offset;
# } Vector;
class c_Vector(Structure):
    _fields_ = [
        ("interface",c_void_p),
        ("p", c_uint),
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
    # Vector * initializeVectorGeneric(uint p, uint64 * vector_container, uint64 * memory, uint dimension, uint offset);
    # Vector * constructVectorGeneric(uint p, uint dimension, uint offset);

    # Vector * initializeVector2(uint p, uint64 * vector_container, uint64 * memory, uint dimension, uint offset);
    # Vector * constructVector2(uint p, uint dimension, uint offset);
    # void freeVector(Vector * v); 

    CSteenrod.constructVectorGeneric.argtypes = [c_uint, c_uint, c_uint]    
    CSteenrod.constructVectorGeneric.restype = POINTER(c_Vector)
    
    CSteenrod.constructVector2.argtypes = [c_uint, c_uint, c_uint]    
    CSteenrod.constructVector2.restype = POINTER(c_Vector)

    CSteenrod.freeVector.argtypes = [POINTER(c_Vector)]
    
    #void assignVector(Vector * target, Vector * source);
    CSteenrod.assignVector.argtypes = [POINTER(c_Vector), POINTER(c_Vector)]
    
    #void setVectorToZero(Vector * target);
    CSteenrod.setVectorToZero.argtypes = [POINTER(c_Vector)]

    #void packVector(Vector * target, uint * source);
    CSteenrod.packVector.argtypes = [POINTER(c_Vector), POINTER(c_uint)]
    
    #void unpackVector(uint * target, Vector * source);
    CSteenrod.unpackVector.argtypes = [POINTER(c_uint), POINTER(c_Vector)]
    
    #uint getVectorEntry(Vector * v, uint index);
    CSteenrod.getVectorEntry.argtypes = [POINTER(c_Vector), c_uint]    
    CSteenrod.getVectorEntry.restype = c_uint
    
    #void setVectorEntry(Vector * v, uint index, uint value);
    CSteenrod.setVectorEntry.argtypes = [POINTER(c_Vector), c_uint, c_uint]    
    
    # void sliceVector(Vector *result, Vector *source, uint start, uint end);
    # VectorIterator getVectorIterator(Vector * v); 
    # VectorIterator stepVectorIterator(VectorIterator);
    CSteenrod.sliceVector.argtypes = [POINTER(c_Vector), POINTER(c_Vector), c_uint, c_uint]

    CSteenrod.getVectorIterator.argtypes = [POINTER(c_Vector)]
    CSteenrod.getVectorIterator.restype = c_VectorIterator
    CSteenrod.stepVectorIterator.argtypes = [c_VectorIterator]
    CSteenrod.stepVectorIterator.restype = c_VectorIterator

    #void addBasisElementToVectorGeneric(Vector * v, uint idx, uint c);
    CSteenrod.addBasisElementToVectorGeneric.argtypes = [POINTER(c_Vector), c_uint, c_uint]  
    
    #void addVectorsGeneric(Vector * target, Vector * source, uint c);
    CSteenrod.addVectorsGeneric.argtypes = [POINTER(c_Vector), POINTER(c_Vector), c_uint]  
    
    #void scaleVectorGeneric(Vector * v, uint c);
    CSteenrod.scaleVectorGeneric.argtypes = [POINTER(c_Vector), c_uint] 

    #void addBasisElementToVector2(Vector * v, uint idx, uint c);
    CSteenrod.addBasisElementToVector2.argtypes = [POINTER(c_Vector), c_uint, c_uint]
    
    #void addVectors2(Vector * target, Vector * source, uint c);
    CSteenrod.addVectors2.argtypes = [POINTER(c_Vector), POINTER(c_Vector), c_uint]  
    
    #void scaleVector2(Vector * v, uint c);
    CSteenrod.scaleVector2.argtypes = [POINTER(c_Vector), c_uint] 
    
    
    #uint vectorToString(char * buffer, Vector * v);
    CSteenrod.vectorToString.argtypes = [c_char_p, POINTER(c_Vector)]
    CSteenrod.vectorToString.restype = c_uint
    
    #Vector* constructMatrixGeneric(uint p, uint rows, uint cols);
    CSteenrod.constructMatrixGeneric.argtypes = [c_uint, c_uint, c_uint]
    CSteenrod.constructMatrixGeneric.restype = POINTER(c_Matrix)
    
    #Vector* constructMatrix2(uint p, uint rows, uint cols);
    CSteenrod.constructMatrix2.argtypes = [c_uint, c_uint, c_uint]
    CSteenrod.constructMatrix2.restype = POINTER(c_Matrix)


    #void rowReduce(Vector **matrix, int * column_to_pivot_row, uint rows);
    CSteenrod.rowReduce.argtypes = [POINTER(c_Matrix), POINTER(c_int)]