import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,os.pardir,"python"))
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,os.pardir,"python", "CWrappers"))

#import pytest
import cAlgebra 

A3 = cAlgebra.getAlgebra("AdemAlgebra", p=3, max_degree=200)

def test_Adem_exhaustive(max_deg):
    import random
    for total in range(max_deg):
        for x_deg in range(total):
            y_deg = total - x_deg
            x_basis = A3.py_algebra.basis(x_deg)
            y_basis = A3.py_algebra.basis(y_deg)
            for x in x_basis:
                for y in y_basis:
                    py_prod = x*y
                    c_prod = A3.multiply(x,y)
                    if(py_prod != c_prod):
                        print("Discrepancy (%s, %s):" % (x, y))
                        print("    py_prod is:", py_prod)
                        print("    c_prod is:", c_prod)
                        print("    difference is:",  (py_prod + (3-1) * c_prod))
                        print("")
                        input("PRESS ENTER TO CONTINUE.")
                    if(random.randint(0,10000)==0):
                        print(c_prod)

test_Adem_exhaustive(80)