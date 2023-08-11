

# this is an example on how to load a localized hamiltonian and compare it's spectrum to the original spectrum
# coeficients h of the localized hamiltonians are stored in the localized folder

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigvalsh
from pauli_strings import local2

def build_local(N, h, taus):
    """compute sum_i h_i tau_i"""
    d = 2**N
    H = np.zeros((d,d), dtype=taus[0].dtype)
    for i in range(len(taus)):
        H[taus[i].row, taus[i].col] += h[i]*taus[i].data
    return H


N = 12
k = 0
taus = local2(N) #list of 2-local pauli strings
h = np.loadtxt('localized/{}/{}'.format(N,k)) # load the kth realization of size N
Hlocal = build_local(N,h,taus) # compute sum_i h_i tau_i
E0 = np.loadtxt('goe_spectra/{}/{}'.format(N,k)) # load the original spectrum
e0 = eigvalsh(Hlocal) # diagonalize the localized hamiltonian

# plot both spectra
plt.hist(E0,bins=100,histtype='step', label='GOE')
plt.hist(e0,bins=100,histtype='step', label='localized')
plt.legend()
plt.show()
