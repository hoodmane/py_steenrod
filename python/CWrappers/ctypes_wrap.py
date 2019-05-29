from ctypes import *
CSteenrod = cdll.LoadLibrary("C/lib/libCSteenrod.so")
from ctypes_combinatorics import *
from ctypes_FpVector import *
from ctypes_algebra import *
from ctypes_milnor import *
from ctypes_modules import *
from ctypes_resolution import *

wrap_combinatorics(CSteenrod)
wrap_FpVector(CSteenrod)
wrap_algebra(CSteenrod)
wrap_modules(CSteenrod)
wrap_milnor(CSteenrod)
wrap_resolution(CSteenrod)
 

