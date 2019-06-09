import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,os.pardir,"python"))
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,os.pardir,"python","CWrappers"))

import json
import pytest
from cResolution import *
external_dims_dict = {}

with open('output_data/S2-100-external.json') as external_json:  
    external_dims_dict[2] = json.load(external_json)

with open('output_data/S3-100-external.json') as external_json:  
    external_dims_dict[3] = json.load(external_json)

def getDimensions(cRes):
    buffer = (c_char * 1000000)()
    CSteenrod.Resolution_gradedDimensionString(buffer, cRes)
    charp = c_char_p(buffer[:])
    s = charp.value.decode("utf-8")
    return json.loads(s) 

@pytest.mark.parametrize("p, degree", [(2,50),(3,75)])
@pytest.mark.parametrize("algebra_type", ["AdemAlgebra","MilnorAlgebra"])
def test_resolution(algebra_type, p, degree):
    A = cAlgebra.getAlgebra(algebra_type, p=p, max_degree=degree)
    M = steenrod_module.FiniteSteenrodModule(p=p)
    x0 = M.add_basis_element("x0", 0)
    M.validate()
    cM = cFiniteDimensionalModule.toC(M, algebra_type)
    res = resolve(M, degree)
    our_dims = getDimensions(res)
    external_dims = external_dims_dict[p]
    for i in range(len(our_dims)):
        for j in range(len(our_dims[i])):
           if external_dims[i][j] != our_dims[i][j]:
               print(i, j, "external_dims(i,j):",external_dims[i][j]," our_dims(i,j):",our_dims[i][j])
               assert False