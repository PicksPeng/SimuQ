import numpy as np

from simuq.environment import qubit
from simuq.hamiltonian import TIHamiltonian
from simuq.qsystem import QSystem

qs = QSystem()
n = 12
ql = [qubit(qs) for i in range(n)]
link = [(i, (i + 1) % n) for i in range(n)]

h0 = TIHamiltonian.empty(n)
for (q0, q1) in link:
    h0 += ql[q0].Z() * ql[q1].Z()


# h1 = TIHamiltonian.empty(n)
# for i in range(n):
#     h1 += ql[i].X()


qs.add_evolution(h0, 1)
# qs.add_evolution(h1, 1)
