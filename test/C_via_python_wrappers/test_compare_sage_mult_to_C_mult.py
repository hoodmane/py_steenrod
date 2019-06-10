import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,os.pardir,"python"))
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,os.pardir,"python", "CWrappers"))

import json
import pytest

import cAlgebra 

p_max_deg_pairs = [(2, 40), (3, 100), (5,200), (7,400)]
sage_products_dict = {}
for algebra_type in ["Adem", "Milnor"]:
    for (p, max_deg) in p_max_deg_pairs:
        with open('output_data/sage_steenrod_mults/sage_steenrod_mults_%s_%s_max_%s.json' % (p, algebra_type, max_deg)) as external_json:  
            sage_products_dict[(algebra_type, p)] = json.load(external_json)

def basis_elt_to_tuples(k):
    if len(k)>0 and type(k[0]) == list:
        k[0] = tuple(k[0])
        k[1] = tuple(k[1])
    return tuple(k)


@pytest.mark.parametrize("algebra_type", ["Milnor", "Adem"])
@pytest.mark.parametrize("p,max_deg", p_max_deg_pairs)
def test_Adem_exhaustive(algebra_type, p, max_deg):
    sage_products = sage_products_dict[(algebra_type, p)]
    A = cAlgebra.getAlgebra(algebra_type + "Algebra", p=p, max_degree=max_deg)
    for degree_d_products in sage_products:
        for entry in degree_d_products:
            if(len(entry[0]) == 0 or len(entry[1])==0):
                continue
            x = A.py_algebra.get_basis_element(basis_elt_to_tuples(entry[0]))
            y = A.py_algebra.get_basis_element(basis_elt_to_tuples(entry[1]))
            res = A.multiply(x,y)
            sage_res_dict = {}
            for k,v in entry[2]:
                k = basis_elt_to_tuples(k)
                sage_res_dict[k] = v
            sage_res = A.py_algebra.get_element(sage_res_dict)
            assert res == sage_res
                         

