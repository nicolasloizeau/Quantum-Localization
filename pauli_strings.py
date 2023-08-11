import os
from scipy.sparse import coo_matrix, kron
import numpy as np
from itertools import product, combinations



paulis = [np.eye(2), np.array([[0,1],[1,0]]), 1j*np.array([[0,-1],[1,0]]), np.array([[1,0],[0,-1]])]
paulis_sparse = [coo_matrix(p, dtype='complex128') for p in paulis]

def operator_from_indexes(indexes, dtype='float64'):
    """
    indexes : list of pauli string indexes (eg [0,1,2,0,3])
    return : coo_matrix representing a pauli string (eg 1XY1Z)
    """
    op = paulis_sparse[indexes[0]]
    for i in indexes[1:]:
        op = kron(op, paulis_sparse[i], format='coo')
    if dtype=='float64':
        op = op.real
    return coo_matrix(op, dtype=dtype)



def local2(N, dtype='float64', sites=[]):
    """
    generates a list of 2-local pauli strings of lenght N
    """
    taus = []
    # 1-local
    for i in range(N):
        pauli_indexes = np.zeros(N, dtype=int)
        pauli_indexes[i] = 3
        taus.append(operator_from_indexes(pauli_indexes, dtype=dtype))
    # 2-local
    for i, j in list(combinations(range(N), 2)):
        for k, l in product(range(1,4), repeat=2):
            # exclude complex strings if dtype=float64
            if dtype=='complex128' or not((k==2 and l!=2) or (l==2 and k!=2)):
                pauli_indexes = np.zeros(N, dtype=int)
                pauli_indexes[i] = k
                pauli_indexes[j] = l
                taus.append(operator_from_indexes(pauli_indexes, dtype=dtype))
    return taus



def local1(N, dtype='float64'):
    """
    generates a list of single qubit Z pauli strings
    (eg Z1111, 1Z111, 11Z11 ....)
    """
    taus = []
    for i in range(N):
        pauli_indexes = np.zeros(N, dtype=int)
        pauli_indexes[i] = 3
        taus.append(operator_from_indexes(pauli_indexes, dtype=dtype))
    return taus
