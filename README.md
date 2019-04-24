### Steenrod Modules

Requires Python 3 (tested in python 3.6.7). 
To use it you should download the repository and then run `./steenrod_modules` in the root directory.

You can load a module from a file by saying `M = SteenrodModule.from_Bruner_file(file_path)`,
for instance `file_path` could be `"Bruner-modules/CP-1"`. You can write a module to a file with `M.to_Bruner_file(file_path)`.

The tensor product of modules `M` and `N` is `M * N`, the dual of a module `M` is `~M`, and the truncation is `M[min:max]`. 
You can create a new empty module with 
`M = SteenrodModule(p = a_prime)`
and then add basis vectors by saying
`x = M.add_basis_element("name", degree)`
and add relations with
`M.add_Sq_action(Sq_degree, basis_element, output_vector)`
You can check the relations define a correct action with `M.validate()`.

So for instance:
```python
A = AdemAlgebra(p=2)
M = SteenrodModule(p=2)
x0 = M.add_basis_element("x0", 0)
x1 = M.add_basis_element("x1", 1)
y1 = M.add_basis_element("y1", 1)
M.add_Sq_action(1, x0, x1 + y1)
M.validate()
```
Now `A.Sq(1) * x0` will return `x1 + y1`.

You could build the same module by saying:
```python
M = SteenrodModule(p=2)
x0 = M.add_basis_element("x0", 0)
x1 = M.add_basis_element("x1", 1)
M.add_Sq_action(1, x0, x1)
N = (M*M)[0:1]
x0x0 = N.get_basis_element("x0*x0")
```
though now the elements will be named `x0*x0`, `x0*x1`, and `x1*x0`. Now `A.Sq(1) * x0x0` will return `x0*x1 + x1*x0`.

`M.get_failed_relations()` will return a list of the Adem relations that are violated in `M`.
