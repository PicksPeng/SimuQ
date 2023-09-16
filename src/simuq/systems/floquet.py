import numpy as np

from simuq.environment import Qubit
from simuq.qsystem import QSystem

# The Floquet Hamiltonian system in arXiv:2107.07311

L = 5
n_repetition = 5
Js = [
    2 * np.pi * 0.01084,
    2 * np.pi * 0.00028,
]  # We assume they are transverse here, which can be specified in detail.
l = len(Js)
geps = 0.1
hs = [0.1, 0.2, 0.3, 0.4, 0.5]
t1 = 0.1
t2 = 0.2
t3 = 0.5

FloquetQS = QSystem()
qs = FloquetQS
ql = [Qubit(qs) for i in range(L)]

hflip = 0
for i in range(L):
    hflip += geps / 2 * ql[i].X

hdis = 0
for i in range(L):
    hdis += hs[i] / 2 * ql[i].Z

hint = 0
for j in range(1, l + 1):
    for i in range(L - j):
        hint += Js[j - 1] / 2 * (ql[i].X * ql[i + j].X + ql[i].Y * ql[i + j].Y)

for i in range(n_repetition):
    qs.add_evolution(hflip, t1)
    qs.add_evolution(hdis, t2)
    qs.add_evolution(hflip, t1)
    qs.add_evolution(hdis, t2)
    qs.add_evolution(hint, t3)
