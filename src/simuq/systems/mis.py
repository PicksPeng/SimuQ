import numpy as np

from simuq.environment import Qubit
from simuq.qsystem import QSystem


def GenQS(topo="Chain", k=3, dis_num=10, dis_partial=None, Delta=1, Omega=4, alpha=4, T=1):
    def adiabatic(h0, h1):
        def f(t):
            return h0 + (-1 + 2.0 * t / T) * h1

        return f

    qs = QSystem()

    if topo == "Chain":
        # Chain
        n = k
        link = [(i, (i + 1) % n) for i in range(n - 1)]
    elif topo == "Cycle":
        # Cycle
        n = k
        link = [(i, (i + 1) % n) for i in range(n)]
    elif topo == "Grid":
        # Grid
        n = k * k
        link = []
        for i in range(k):
            for j in range(k):
                if j < k - 1:
                    link.append((i * k + j, i * k + j + 1))
                if i < k - 1:
                    link.append((i * k + j, (i + 1) * k + j))

    ql = [Qubit(qs) for i in range(n)]
    noper = [(ql[i].I - ql[i].Z) / 2 for i in range(n)]
    Omega = 4
    alpha = 4
    h0 = 0
    for q0, q1 in link:
        h0 += alpha * noper[q0] * noper[q1]
    for i in range(n):
        h0 += Omega / 2 * ql[i].X

    Delta = 1
    h1 = 0
    for i in range(n):
        h1 += -Delta * noper[i]

    if dis_partial is None:
        dis_partial = dis_num
    qs.add_time_dependent_evolution(
        adiabatic(h0, h1), np.linspace(0, T, dis_num + 1)[: dis_partial + 1]
    )

    return qs
