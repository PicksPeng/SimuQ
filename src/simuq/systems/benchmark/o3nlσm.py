# Generate an 1+1D O(3) NLσM model
# The following is described in Phys. Rev. A 107, 042404
# H = Jx Σ_{(j, k)} σ_{j, k}·σ_{(j+1), k} + Jy Σ_{(j, k)} σ_{j, k}·σ_{j, (k+1)}
# Here σ_{j, k}=(X_{j, k}, Y_{j, k}, Z_{j, k})

from simuq.environment import Qubit
from simuq.qsystem import QSystem


def GenQS(n, m, T=1, Jx=0.5, Jy=0.8):
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
