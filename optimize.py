
# this is the main localization code,
# followed by an example of localizing a 8 qubit GOE spectrum

import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from numpy.linalg import norm
from scipy.linalg import eigh, eigvalsh


class Localizer:

    def __init__(self, E0, taus, maxiter=100, gtol=1e-6, h0=[], disp=True):
        self.N = int(np.log2(len(E0)))
        self.taus = taus
        self.iter = 0
        self.history = []
        self.maxiter = maxiter
        self.gtol = gtol
        self.E0 = E0
        self.h0 = h0
        if len(self.h0) == 0:
            self.h0 = (np.random.random(len(self.taus))-0.5)/len(self.taus)
        self.disp = disp
        self.current_cost = 1


    def build_local(self, h):
        """compute sum_i h_i tau_i"""
        d = 2**self.N
        H = np.zeros((d,d), dtype=self.taus[0].dtype)
        for i in range(len(self.taus)):
            H[self.taus[i].row, self.taus[i].col] += h[i]*self.taus[i].data
        return H

    def cost(self, h):
        E = eigvalsh(self.build_local(h))
        c = norm(E-self.E0)**2
        self.current_cost = c
        return c


    def callback(self, h):
        c = self.current_cost
        self.history.append(c)
        self.iter += 1
        if self.disp:
            print(self.iter, c)

    def localize(self):
        res = minimize(self.cost, self.h0, method='BFGS',
        callback=self.callback, jac=self.gradient,
        options={'gtol':1e-16, 'maxiter':self.maxiter, 'disp':True}, tol=1e-16)
        return res.x


    def gradient(self, h):
        e,v = eigh(self.build_local(h))
        H2 = v@np.diag(self.E0)@v.T
        grad = np.zeros(len(h))
        for i in range(len(h)):
            grad[i] = 2*(h[i]*2**self.N-np.sum(H2[self.taus[i].row, self.taus[i].col]*self.taus[i].data))
        return grad




if __name__ == "__main__":
    from pauli_strings import*
    # this is an example of localizing a 8 qubit GOE spectrum
    N = 8 # system size (qbit)
    k = 0 # file index (there are 200 spectra for each N stored in goe_spectra, here we open the first one)
    iter = 500 # number of iterations
    taus = local2(N) # list of 2-local pauli strings
    E0 = np.loadtxt('goe_spectra/{}/{}'.format(N,k)) # open the spectrum
    # localize the spectrum. The coeficients of the 2-localized H are stored in h
    loc = Localizer(E0, taus, maxiter=iter)
    h = loc.localize()
    # cost vs iteration plot
    plt.loglog(loc.history)
    plt.xlabel('iteration')
    plt.ylabel('cost')
    plt.show()
