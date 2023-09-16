# Generate a transverse field cycle Ising model
# H = J \sum_{j=1}^{n} Z_jZ_{j+1} + h \sum_{j=1}^n X_j  (assuming Z_{n+1}=Z_1)

from simuq.environment import Qubit
from simuq.qsystem import QSystem


def GenQS(n, T, J, h):
    qs = QSystem()
    q = [Qubit(qs) for i in range(n)]
    H = 0
    for i in range(n):
        H = H + J * q[i].Z * q[(i + 1) % n].Z
    for i in range(n):
        H = H + h * q[i].X
    qs.add_evolution(H, T)
    return qs
