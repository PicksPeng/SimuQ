import numpy as np
from simuq.environment import qubit
from simuq.hamiltonian import TIHamiltonian
from simuq.qsystem import QSystem

T = 30

def adiabatic(h0, h1) :
    def f(t) :
        return h0 + (-1 + 2. * t / T) * h1
    return f
m = 5

qs = QSystem()
n = 5
ql = [qubit(qs) for i in range(n)]
link = [(i, (i + 1) % n) for i in range(n-1)]
Omega = 1
h0 = TIHamiltonian.empty(n)
for (q0, q1) in link :
    h0 += (ql[q0].I() - ql[q0].Z()) * (ql[q1].I() - ql[q1].Z())
    #h += ql[q0].Z() * ql[q1].Z()
for i in range(n) :
    h0 += Omega / 2 * ql[i].X()

Delta = 5
h1 = TIHamiltonian.empty(n)
for i in range(n) :
    h1 += - Delta * (ql[i].I() - ql[i].Z())


qs.add_time_dependent_evolution(adiabatic(h0, h1), np.linspace(0, T, m))
