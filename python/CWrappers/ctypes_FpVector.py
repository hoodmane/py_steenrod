from ctypes import *

class c_Vector(Structure):
    _fields_ = []


def wrap_FpVector(CSteenrod):        
    #Vector * constructVector(uint p, uint dimension);
    CSteenrod.constructVector.argtypes = [c_uint, c_uint]    
    CSteenrod.constructVector.restype = POINTER(c_Vector)
    
    #void freeVector(Vector * v);
    CSteenrod.freeVector.argtypes = [POINTER(c_Vector)]
    
    #void assignVector(uint p, Vector * target, Vector * source);
    CSteenrod.assignVector.argtypes = [c_uint, POINTER(c_Vector), POINTER(c_Vector)]
    
    #void packVector(uint p, Vector * target, uint * source);
    CSteenrod.packVector.argtypes = [c_uint, POINTER(c_Vector), POINTER(c_uint)]
    
    #void unpackVector(uint p, uint * target, Vector * source);
    CSteenrod.unpackVector.argtypes = [c_uint, POINTER(c_uint), POINTER(c_Vector)]
    
    #uint vectorToString(char * buffer, uint p, Vector * v);
    CSteenrod.vectorToString.argtypes = [c_char_p, c_uint, POINTER(c_Vector)]
    CSteenrod.vectorToString.restype = c_uint
    
    #uint getVectorEntry(uint p, Vector * v, uint index);
    CSteenrod.getVectorEntry.argtypes = [c_uint, POINTER(c_Vector), c_uint]    
    CSteenrod.getVectorEntry.restypes = c_uint
    
    #void setVectorEntry(uint p, Vector * v, uint index, uint value);
    CSteenrod.setVectorEntry.argtypes = [c_uint, POINTER(c_Vector), c_uint, c_uint]    
    
    #void addBasisElementToVectorGeneric(uint p, Vector * v, uint idx, uint c);
    CSteenrod.addBasisElementToVectorGeneric.argtypes = [c_uint, POINTER(c_Vector), c_uint, c_uint]  
    
    #void addVectorsGeneric(uint p, Vector * target, Vector * source, uint c);
    CSteenrod.addVectorsGeneric.argtypes = [c_uint, POINTER(c_Vector), POINTER(c_Vector), c_uint]  
    
    #void scaleVectorGeneric(uint p, Vector * v, uint c);
    CSteenrod.scaleVectorGeneric.argtypes = [c_uint, POINTER(c_Vector), c_uint] 

    #void addBasisElementToVector2(uint p, Vector * v, uint idx, uint c);
    CSteenrod.addBasisElementToVector2.argtypes = [c_uint, POINTER(c_Vector), c_uint, c_uint]
    
    #void addVectors2(uint p, Vector * target, Vector * source, uint c);
    CSteenrod.addVectors2.argtypes = [c_uint, POINTER(c_Vector), POINTER(c_Vector), c_uint]  
    
    #void scaleVector2(uint p, Vector * v, uint c);
    CSteenrod.scaleVector2.argtypes = [c_uint, POINTER(c_Vector), c_uint] 
    
    
    
    
    
