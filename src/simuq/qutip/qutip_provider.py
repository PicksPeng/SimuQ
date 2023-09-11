import numpy as np
import qutip as qp

from simuq.provider import BaseProvider


class QuTiPProvider(BaseProvider):
    def __init__(self):
        self.prog = None
        self.fin = None

    def compile(self, qs, initial_state=None, verbose=0):
        self.n = qs.num_sites
        if initial_state is None:
            self.init = qp.basis(1 << self.n)
        else:
            self.init = initial_state
        self.prog = (qs.to_qutip(), qs.total_time())
        self.qs_names = qs.print_sites()
        if verbose >= 0:
            print("Compiled.")
        # return self.prog

    def evaluate_Hamiltonian(self, t):
        if self.prog is None:
            raise Exception("No compiled job in record.")
        M = 0
        for i in range(len(self.prog[0])):
            (H, f) = self.prog[0][i]
            M += H * f(t, None)
        return M

    def run(self, shots=None, on_simulator=None, verbose=0):
        if self.prog is None:
            raise Exception("No compiled job in record.")
        self.fin = qp.sesolve(self.prog[0], self.init, [0, self.prog[1]])
        if verbose >= 0:
            print("Solved.")
        # return self.fin

    def results(self, verbose=0):
        if self.fin is None:
            raise Exception("No submitted job in record.")
        self.res = dict()
        for i in range(1 << self.n):
            self.res[QuTiPProvider.to_bin(i, self.n)] = np.abs(self.fin.states[1][i][0][0]) ** 2
        return self.res
