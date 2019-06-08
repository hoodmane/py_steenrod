import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,os.pardir,"python"))
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,os.pardir,"python", "CWrappers"))

#import pytest
print("ho")
import cMilnorAlgebra 

A3 = cMilnorAlgebra.construct(p=3, degree=200)

def test_milnor_exhaustive(max_deg):
    import random
    for x_deg in range(max_deg):
        for y_deg in range(max_deg):
            x_basis = cMilnorAlgebra.getBasis(A3, x_deg)
            y_basis = cMilnorAlgebra.getBasis(A3, y_deg)
            for x in x_basis:
                for y in y_basis:
                    py_prod = x*y
                    c_prod = cMilnorAlgebra.multiply(x,y)
                    if(py_prod != c_prod):
                        print("Discrepancy (%s, %s):" % (x, y))
                        print("    py_prod is:", py_prod)
                        print("    c_prod is:", c_prod)
                        print("    difference is:",  (py_prod + (3-1) * c_prod))
                        print("")
                        input("PRESS ENTER TO CONTINUE.")
                    if(random.randint(0,10000)==0):
                        print(c_prod)

#test_milnor_exhaustive()