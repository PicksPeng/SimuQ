import numpy as np
from simuq.qsystem import QSystem
from simuq.environment import qubit
from simuq.hamiltonian import TIHamiltonian

qs = QSystem()
n = 3
ql = [qubit(qs) for i in range(n)]
link = [(0, 1), (1, 2)]
T = 1
h = TIHamiltonian.empty(n)
for (q0, q1) in link :
    h += (ql[q0].I() - ql[q0].Z()) * (ql[q1].I() - ql[q1].Z())
qs.add_evolution(h, T)
