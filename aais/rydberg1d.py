from simuq.environment import qubit, fock
from simuq.qmachine import *
from simuq.expression import Expression
from simuq.hamiltonian import TIHamiltonian
import numpy as np

Rydberg = QMachine()

C6 = 862690 * 2 * np.pi
#C6 = 4

n = 5

q = [qubit(Rydberg) for i in range(n)]

l = (C6 / 4) ** (1. / 6)
x = [0] + [GlobalVar(Rydberg, init_value = l * i) for i in range(1, n)]

sys_h = TIHamiltonian.empty(n)
for i in range(n) :
    for j in range(i) :
        sys_h += (C6 / 4 / (x[i] - x[j])**6) * (q[j].I() - q[j].Z()) * (q[i].I() - q[i].Z())
Rydberg.set_sys_ham(sys_h)

for i in range(n) :
    L = SignalLine(Rydberg)
    ins = Instruction(L, 'native', f'Detuning of site {i}')
    d = LocalVar(ins)
    ins.set_ham(- d * (q[i].I() - q[i].Z()))

for i in range(n) :
    L = SignalLine(Rydberg)
    ins = Instruction(L, 'native', f'Rabi of site {i}')
    o = LocalVar(ins)
    p = LocalVar(ins)
    ins.set_ham(o / 2 * (Expression.cos(p) * q[i].X() - Expression.sin(p) * q[i].Y()))