# Generate an electron-photon dynamics model
# The following is described in New J. Phys. 24, 9, 093017
# H = Hel + Hph + Hep
# Hel = 1/2 Σ_j eps_j Z_j + 1/2 Σ_k Σ_j V_jk (X_j X_k + Y_j Y_k)
# Hph =

from simuq.environment import Qubit
from simuq.qsystem import QSystem


def GenQS(n, m, Jx, Jy):
    qs = QSystem()
    q = [[Qubit(qs) for j in range(m)] for i in range(n)]
    H = 0
    for i in range(n - 1):
        for j in range(m):
            H += Jx * q[i][j].X * q[i + 1][j].X
            H += Jx * q[i][j].Y * q[i + 1][j].Y
            H += Jx * q[i][j].Z * q[i + 1][j].Z
    for i in range(n):
        for j in range(m - 1):
            H += Jx * q[i][j].X * q[i][j + 1].X
            H += Jx * q[i][j].Y * q[i][j + 1].Y
            H += Jx * q[i][j].Z * q[i][j + 1].Z
    qs.add_evolution(H, T)
    return qs
