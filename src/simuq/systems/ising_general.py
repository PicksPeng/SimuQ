import numpy as np
from scipy.constants import h

from simuq.environment import Qubit
from simuq.qsystem import QSystem

n = 3
m = 3
T = 1
tmp = np.random.uniform(0, 1, (n, n))
J = (tmp + tmp.T) / 2
qs = QSystem()
qubits = [Qubit(qs) for i in range(n)]
H0 = 0
H1 = 0
for i in range(n):
    for j in range(n):
        H0 += -J[i][j] * qubits[i].Z * qubits[j].Z

for i in range(n):
    H0 += -h * qubits[i].Z
    H1 += qubits[i].X


def ising_model(H0, H1, Tau, T):
    def f(t):
        return H0 - Tau(t / T) * H1

    return f


qs.add_time_dependent_evolution(ising_model(H0, H1, lambda t: 1 - t**2, T), np.linspace(0, T, m))
