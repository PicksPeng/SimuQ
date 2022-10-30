import numpy as np

from simuq.environment import qubit
from simuq.hamiltonian import Empty
from simuq.qsystem import QSystem

T = 1


def adiabatic(h0, h1):
    def f(t):
        return h0 + (-1 + 2.0 * t / T) * h1

    return f

qs = QSystem()
n = 3
ql = [qubit(qs) for i in range(n)]
link = [(i, (i + 1) % n) for i in range(n - 1)]
noper = [(ql[i].I - ql[i].Z) / 2 for i in range(n)]
Omega = 5
alpha = 4
h0 = Empty
for (q0, q1) in link:
    h0 += alpha * noper[q0] * noper[q1]
    # h += ql[q0].Z * ql[q1].Z
for i in range(n):
    h0 += Omega / 2 * ql[i].X

Delta = 4
h1 = Empty
for i in range(n):
    h1 += -Delta * noper[i]

def set_qs(m = 10, step = 10) :
    qs.clear_evos()
    qs.add_time_dependent_evolution(adiabatic(h0, h1), np.linspace(0, T, m + 1)[:step + 1])

set_qs()
