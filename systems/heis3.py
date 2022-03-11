import numpy as np
from simuq.qsystem import QSystem
from simuq.environment import qubit
from simuq.hamiltonian import TIHamiltonian

qs = QSystem()
n = 3
ql = [qubit(qs) for i in range(n)]
link = [(0, 1), (1, 2)]
T = np.pi / 4
h = TIHamiltonian.empty(n)
for (q0, q1) in link :
    h += ql[q0].X() * ql[q1].X()
    h += ql[q0].Y() * ql[q1].Y()
    h += ql[q0].Z() * ql[q1].Z()
qs.add_evolution(h, T)
