# Generate a Kitaev chain coupled to a Z2 gauge field
# The following comes from Aramthottil et al., Phys. Rev. B, L041101, 2022
# H = μ/2 Σ_{j=1}^{n-1} Z_jZ_{j+1} - Σ_{j=1}^n (t X_j + h Z_j)

from simuq.environment import Qubit
from simuq.qsystem import QSystem


def GenQS(n, T=1, mu=1, t=0.3, h=0.5):
    qs = QSystem()
    q = [Qubit(qs) for i in range(n)]
    H = 0
    for i in range(n - 1):
        H += mu / 2 * q[i].Z * q[i + 1].Z
    for i in range(n):
        H += -t * q[i].X - h * q[i].Z
    qs.add_evolution(H, T)
    return qs
