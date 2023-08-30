import numpy as np

from simuq.environment import qubit
from simuq.qsystem import QSystem

qs = QSystem()
n = 3
ql = [qubit(qs) for i in range(n)]
link = [(i, i + 1) for i in range(n - 1)]
T = 1.0
scaler = 1000
h = 0
for i in range(n):
    h += scaler * ql[i].X
qs.add_evolution(h, T / scaler)
