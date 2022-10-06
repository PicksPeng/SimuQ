from simuq.environment import qubit, fock
from simuq.qmachine import *
from simuq.expression import Expression
from simuq.hamiltonian import Empty
import numpy as np

Rydberg = QMachine()

C6 = 862690 * 2. * np.pi
#C6 = 4.

n = 3

q = [qubit(Rydberg) for i in range(n)]

l = (C6 / 4) ** (1. / 6) / (2 - 2 * np.cos(2 * np.pi / n)) ** 0.5
x = [(0, 0.)] + [(GlobalVar(Rydberg, init_value = l * (np.cos(i * 2 * np.pi / n) - 1)), GlobalVar(Rydberg, init_value = l * np.sin(i * 2 * np.pi / n))) for i in range(1, n)]
noper = [(q[i].I - q[i].Z) / 2 for i in range(n)]

sys_h = Empty
for i in range(n) :
    for j in range(i) :
        dsqr = (x[i][0] - x[j][0])**2 + (x[i][1] - x[j][1])**2
        sys_h += (C6 / 4 / (dsqr ** 3)) * noper[i] * noper[j]
Rydberg.set_sys_ham(sys_h)

for i in range(n) :
    L = SignalLine(Rydberg)
    ins = Instruction(L, 'native', f'Detuning of site {i}')
    d = LocalVar(ins)
    ins.set_ham(- d * noper[i])

for i in range(n) :
    L = SignalLine(Rydberg)
    ins = Instruction(L, 'native', f'Rabi of site {i}')
    o = LocalVar(ins)
    p = LocalVar(ins)
    ins.set_ham(o / 2 * (Expression.cos(p) * q[i].X - Expression.sin(p) * q[i].Y))
