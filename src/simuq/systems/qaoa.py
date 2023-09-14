import numpy as np

from simuq.environment import Qubit
from simuq.qsystem import QSystem


def GenQS(n=12, p=3):
    qs = QSystem()
    ql = [Qubit(qs) for i in range(n)]
    link = [(i, (i + 1) % n) for i in range(n)]
    parameter_list = (
        np.array(
            [
                0.5702193 * 2,
                -0.58631086,
                0.85160685 * 2,
                -1.7058538,
                0.29468536 * 2,
                -1.132814,
            ]
        )
        / 2
    )
    """
    parameter_list = np.random.uniform(0, 1, 2*p)
    """
    for i in range(p):
        h = 0
        for q0, q1 in link:
            h += parameter_list[2 * i] * ql[q0].Z * ql[q1].Z
        qs.add_evolution(h, 1)
        h = 0
        for q0 in range(n):
            h += parameter_list[2 * i + 1] * ql[q0].X
        qs.add_evolution(h, 1)
    return qs
