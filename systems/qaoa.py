import numpy as np

from simuq.environment import qubit
from simuq.hamiltonian import TIHamiltonian
from simuq.qsystem import QSystem

qs = QSystem()
n = 12
ql = [qubit(qs) for i in range(n)]
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
for i in range(3):
    h0 = TIHamiltonian.empty(n)
    for (q0, q1) in link:
        h0 += parameter_list[2 * i] * ql[q0].Z() * ql[q1].Z()
    qs.add_evolution(h0, 1)
    h0 = TIHamiltonian.empty(n)
    for q0 in range(n):
        h0 += parameter_list[2 * i + 1] * ql[q0].X()
    qs.add_evolution(h0, 1)


# h1 = TIHamiltonian.empty(n)
# for i in range(n):
#     h1 += ql[i].X()


# qs.add_evolution(h1, 1)
