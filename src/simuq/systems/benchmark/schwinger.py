# Generate a Schwinger model in pure spin Hamiltonian
# The following is described in [Hamer et al., Phys. Rev. D 56, 55-67, 1997]
# H = m/2 Σ_j (-1)^j Z_n + w/2 Σ_{j=1}^{N-1} (X_jX_{j+1}+Y_jY_{j+1})
#   + J Σ_{j=1}^{N-1} (eps_0 + 1/2 Σ_{k=1}^j (Z_k+(-1)^k))^2

from simuq.environment import Qubit
from simuq.qsystem import QSystem


def GenQS(n, T=1, m=0.5, w=1, J=1, eps_0=0):
    qs = QSystem()
    q = [Qubit(qs) for i in range(n)]
    H = 0
    for j in range(n):
        H += m / 2 * (-1) ** j * q[j].Z
    for j in range(n - 1):
        H += w / 2 * (q[j].X * q[j + 1].X + q[j].Y * q[j + 1].Y)
    for j in range(n - 1):
        S = 0
        for k in range(j):
            S += q[k].Z + (-1) ** k
        H = J * (eps_0 + S) * (eps_0 + S)
    qs.add_evolution(H, T)
    return qs
