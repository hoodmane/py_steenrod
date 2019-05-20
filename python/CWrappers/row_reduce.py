from ctypes import *
from ctypes_wrap import *

    
if __name__ == "__main__":
    p = 7
    CSteenrod.initializePrime(p)
    state = c_row_reduce_state();
    state.p = p
    state.row_capacity = 100
    state.column_capacity = 100
    state.pivot = -1
    state.source_dim = 5
    state.target_dim = 8
    state.matrix = CSteenrod.allocate_matrix(state.row_capacity, state.column_capacity)
    import random
    for row in range(state.source_dim):
        for col in range(state.target_dim):
            state.matrix[row][col] = random.randint(0,p-1)
        state.matrix[row][state.target_dim + row] = 1
    print([r[0:state.source_dim + state.target_dim] for r in state.matrix[0:state.source_dim]])
    
    
    

