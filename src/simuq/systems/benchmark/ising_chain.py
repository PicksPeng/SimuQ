# Generate a transverse field chain Ising model
# H = J Σ_{j=1}^{n-1} Z_jZ_{j+1} + h Σ_{j=1}^n X_j

from simuq.environment import Qubit
from simuq.qsystem import QSystem


def GenQS(n, T, J, h):
    qs = QSystem()
    q = [Qubit(qs) for i in range(n)]
    H = 0
    for i in range(n - 1):
        H = H + J * q[i].Z * q[i + 1].Z
    for i in range(n):
        H = H + h * q[i].X
    qs.add_evolution(H, T)
    return qs
