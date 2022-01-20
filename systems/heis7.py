import numpy as np
from simuq.qsystem import QSystem
from simuq.environment import qubit
from simuq.hamiltonian import TIHamiltonian

qs = QSystem()
n = 7
ql = [qubit(qs) for i in range(7)]
link = [(0, 1), (1, 2), (1, 3), (3, 5), (4, 5), (5, 6)]
T = np.pi
h = TIHamiltonian.empty(n)
for (q0, q1) in link :
    h += ql[q0].X() * ql[q1].X()
    h += ql[q0].Y() * ql[q1].Y()
    h += ql[q0].Z() * ql[q1].Z()
qs.add_evolution(h, np.pi)
