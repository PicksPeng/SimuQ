import numpy as np

from simuq.environment import Qubit
from simuq.qsystem import QSystem

qs = QSystem()
n = 7
ql = [Qubit(qs) for i in range(7)]
link = [(0, 1), (1, 2), (1, 3), (3, 5), (4, 5), (5, 6)]
T = np.pi
h = 0
Jx = 1
Jy = 1
Jz = 1
for q0, q1 in link:
    h += Jx * ql[q0].X * ql[q1].X
    h += Jy * ql[q0].Y * ql[q1].Y
    h += Jz * ql[q0].Z * ql[q1].Z
a = 1
for i in range(n):
    h += a * ql[i].Z
qs.add_evolution(h, T)
