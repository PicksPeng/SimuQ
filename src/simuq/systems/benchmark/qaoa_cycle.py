# The QAOA algorithm with random parameters.
# H1 = Σ_{(j, k)\in E} Z_jZ_k
# H2 = Σ_j X_j
# The algorithm evolves alternativey under H1 and H2,
# whose length follows the parameter_list

import numpy as np

from simuq.environment import Qubit
from simuq.qsystem import QSystem


def GenQS(n=12, p=3, parameter_list=None):
    qs = QSystem()
    q = [Qubit(qs) for i in range(n)]
    link = [(i, (i + 1) % n) for i in range(n)]
    if parameter_list is None:
        parameter_list = np.random.uniform(0, 1, 2 * p)
    for i in range(p):
        h = 0
        for j, k in link:
            h += parameter_list[2 * i] * q[j].Z * q[k].Z
        qs.add_evolution(h, 1)
        h = 0
        for j in range(n):
            h += parameter_list[2 * i + 1] * q[j].X
        qs.add_evolution(h, 1)
    return qs
