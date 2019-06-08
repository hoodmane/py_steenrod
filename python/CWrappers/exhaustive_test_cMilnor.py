import pytest
import cMilnorAlgebra 

A3 = cMilnorAlgebra.construct(p=3, degree=200)

for x_deg in range(85):
    for y_deg in range(85):
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