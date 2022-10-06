import numpy as np
from simuq.environment import qubit
from simuq.hamiltonian import TIHamiltonian
from simuq.qsystem import QSystem

T = 20

def adiabatic(h0, h1) :
    def f(t) :
        return h0 + (-1 + 2. * t / T) * h1
    return f
m = 3

qs = QSystem()
n = 3
ql = [qubit(qs) for i in range(n)]
link = [(i, (i + 1) % n) for i in range(n-1)]
noper = [(ql[i].I - ql[i].Z) / 2 for i in range(n)]
Omega = 1
h0 = TIHamiltonian.empty()
for (q0, q1) in link :
    h0 += noper[q0] * noper[q1]
    #h += ql[q0].Z * ql[q1].Z
for i in range(n) :
    h0 += Omega / 2 * ql[i].X

Delta = 3
h1 = TIHamiltonian.empty()
for i in range(n) :
    h1 += - Delta * noper[i]


qs.add_time_dependent_evolution(adiabatic(h0, h1), np.linspace(0, T, m))
