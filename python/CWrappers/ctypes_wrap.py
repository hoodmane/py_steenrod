from ctypes import *
CSteenrod = cdll.LoadLibrary("C/lib/libCSteenrod.so")
from ctypes_combinatorics import *
from ctypes_algebra import *
from ctypes_milnor_datatypes import *
from ctypes_milnor import *


wrap_combinatorics(CSteenrod)
wrap_algebra(CSteenrod)
wrap_milnor_datatypes(CSteenrod)
wrap_milnor(CSteenrod)
    

