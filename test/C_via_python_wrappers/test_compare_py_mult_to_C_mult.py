import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,os.pardir,"python"))
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,os.pardir,"python", "CWrappers"))

import pytest
import cAlgebra 

@pytest.mark.parametrize("algebra_type", ["MilnorAlgebra", "AdemAlgebra"])
@pytest.mark.parametrize("p,max_deg", [(2, 30), (3, 50)])
def test_Adem_exhaustive(algebra_type, p, max_deg):
    A = cAlgebra.getAlgebra(algebra_type, p=p, max_degree=max_deg)
    for total in range(max_deg):
        for x_deg in range(total):
            y_deg = total - x_deg
            x_basis = A.py_algebra.basis(x_deg)
            y_basis = A.py_algebra.basis(y_deg)
            for x in x_basis:
                for y in y_basis:
                    py_prod = x*y
                    c_prod = A.multiply(x,y)
                    if(py_prod != c_prod):
                        print("Discrepancy (%s, %s):" % (x, y))
                        print("    py_prod is:", py_prod)
                        print("    c_prod is:", c_prod)
                        print("    difference is:",  (py_prod + (3-1) * c_prod))
                        print("")
                        assert False