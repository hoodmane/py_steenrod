from ctypes import *


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
    
