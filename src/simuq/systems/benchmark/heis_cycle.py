# Generate a transverse field cycle Ising model
# H = J \sum_{j=1}^{n} σ_j·σ_{j+1} + h \sum_{j=1}^n X_j  (assuming σ_{n+1}=σ_1)
# Here σ_j=(X_j, Y_j, Z_j)

from simuq.environment import Qubit
from simuq.qsystem import QSystem


def GenQS(n, T, J, h):
    qs = QSystem()
    q = [Qubit(qs) for i in range(n)]
    H = 0
    for i in range(n):
        H = H + J * q[i].X * q[(i + 1) % n].X
        H = H + J * q[i].Y * q[(i + 1) % n].Y
        H = H + J * q[i].Z * q[(i + 1) % n].Z
    for i in range(n):
        H = H + h * q[i].X
    qs.add_evolution(H, T)
    return qs
