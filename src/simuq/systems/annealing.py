import numpy as np

from simuq.environment import Qubit
from simuq.qsystem import QSystem


def anneal(h0, h1, T):
    def f(t):
        return (1 - t / T) * h0 + t / T * h1

    return f


n = 3  # num of qubits
m = 3  # discretization
T = 1  # evolution time

qs = QSystem()
q = [Qubit(qs) for i in range(n)]
h0, h1 = 0, 0
for i in range(n):
    h0 += q[i].X
for i in range(n - 1):
    h1 += q[i].Z * q[i + 1].Z

qs.add_td_evolution(anneal(h0, h1, T), np.linspace(0, T, m))
