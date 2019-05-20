from ctypes import *

#typedef struct {
#    unsigned long p;
#    long source_dim;
#    long target_dim;
#    long last_pivot;
#    long row_capacity;
#    long column_capacity;
#    long ** matrix;
#    bool found_cokernel;
#} row_reduce_state;

class c_row_reduce_state(Structure):
    _fields_ = [
        ("p", c_ulong),
        ("source_dim", c_long),
        ("target_dim", c_long),
        ("pivot", c_long),
        ("row_capacity", c_long),
        ("column_capacity",c_long), 
        ("matrix", POINTER(POINTER(c_long))),
        ("found_cokernel", c_bool),
    ]

def wrap_combinatorics(CSteenrod):
    CSteenrod.initializePrime.argtypes = [c_ulong]
    
    CSteenrod.inverse.argtypes = [c_ulong, c_long]
    CSteenrod.inverse.restype = c_long
    
    CSteenrod.Binomial.argtypes = [c_ulong, c_ulong, c_ulong]
    CSteenrod.Binomial.restype = c_ulong
    CSteenrod.Multinomial.argtypes = [c_ulong, c_ulong, POINTER(c_ulong)]
    CSteenrod.Multinomial.restype = c_ulong
    
    CSteenrod.getXiDegrees.argtypes = [c_ulong]
    CSteenrod.getXiDegrees.restype = POINTER(c_ulong)
    CSteenrod.getTauDegrees.argtypes = [c_ulong]
    CSteenrod.getTauDegrees.restype = POINTER(c_ulong)
    
    CSteenrod.allocate_matrix.argtypes = [c_ulong, c_ulong]
    CSteenrod.allocate_matrix.restype = POINTER(POINTER(c_long))
    CSteenrod.row_reduce.argtypes = [c_row_reduce_state]
