import numpy as np
from simuq.qsystem import QSystem
from simuq.environment import qubit
from simuq.hamiltonian import Empty

def anneal01(h0, h1) :
    def f(t) :
        return t * h0 + (1 - t) * h1
    return f

n = 3 # num of qubits
m = 3 # discretization

qs = QSystem()

ql = [qubit(qs) for i in range(n)]

h0 = Empty
for i in range(n) :
    h0 += ql[i].X

h1 = Empty
for i in range(n - 1) :
    h1 += ql[i].Z * ql[i + 1].Z

qs.add_time_dependent_evolution(anneal01(h0, h1), np.linspace(0, 1, m))
