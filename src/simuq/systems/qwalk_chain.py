from simuq.environment import Qubit
from simuq.qsystem import QSystem

# Quantum walk on a 1D chain
# Antiferromagnetic embedding n = N - 1
# Hpen = \sum_{j=1}^{n-1} Z_j Z_{j+1} + Z_1 + (-1)^n Z_n
# Q = \sum_{j=1}^n X_j


def GenQS(N, lamb, T, encoding="antiferromagnetic"):
    if encoding == "antiferromagnetic":
        n = N - 1
        qs = QSystem()
        ql = [Qubit(qs) for i in range(n)]

        Hpen = 0
        for j in range(n - 1):
            Hpen += ql[j].Z * ql[j + 1].Z
        Hpen += ql[0].Z
        if n % 2 == 0:
            Hpen += ql[n - 1].Z
        else:
            Hpen -= ql[n - 1].Z

        Q = 0
        for j in range(n):
            Q += ql[j].X

        h = lamb * Hpen + Q

        qs.add_evolution(h, T)
        return qs
