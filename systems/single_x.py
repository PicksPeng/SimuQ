import numpy as np
from simuq.qsystem import QSystem
from simuq.environment import qubit
from simuq.hamiltonian import TIHamiltonian

qs = QSystem()
n = 7
ql = [qubit(qs) for i in range(7)]
link = [(0, 1), (1, 2), (1, 3), (3, 5), (4, 5), (5, 6)]
T = 2 * np.pi
h = TIHamiltonian.empty(n)
h += ql[0].X()
qs.add_evolution(h, T)
