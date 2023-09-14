from simuq.environment import Qubit
from simuq.qsystem import QSystem


# Chain or cycle transverse Ising model
def GenQS(n, T, J, h, is_chain=True):
    qs = QSystem()
    q = [Qubit(qs) for i in range(n)]
    H = 0
    for i in range(n - 1 if is_chain else n):
        H = H + J * q[i].Z * q[(i + 1) % n].Z
    for i in range(n):
        H = H + h * q[i].X
    qs.add_evolution(H, T)
    return qs
