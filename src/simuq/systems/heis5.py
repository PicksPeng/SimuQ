import numpy as np

from simuq.environment import Qubit
from simuq.qsystem import QSystem

qs = QSystem()
n = 5
ql = [Qubit(qs) for i in range(n)]
link = [(0, 1), (1, 2), (2, 3), (3, 4)]
T = np.pi / 4
h = 0
for q0, q1 in link:
    h += ql[q0].X * ql[q1].X
    h += ql[q0].Y * ql[q1].Y
    h += ql[q0].Z * ql[q1].Z
for i in range(n):
    h += ql[i].X
qs.add_evolution(h, T)
