import numpy as np

from simuq.environment import fock, qubit
from simuq.expression import Expression
from simuq.hamiltonian import Empty
from simuq.qmachine import *

Rydberg = QMachine()

C6 = 862690 * 2 * np.pi
# C6 = 4

n = 3

q = [qubit(Rydberg) for i in range(n)]

l = (C6 / 4) ** (1.0 / 6)
x = [0] + [GlobalVar(Rydberg, init_value=l * i) for i in range(1, n)]
noper = [(q[i].I - q[i].Z) / 2 for i in range(n)]

sys_h = Empty
for i in range(n):
    for j in range(i):
        sys_h += (C6 / 4 / (x[i] - x[j]) ** 6) * noper[i] * noper[j]
Rydberg.set_sys_ham(sys_h)

# for i in range(n):
L = SignalLine(Rydberg)
ins = Instruction(L, "native", f"Detuning")
d = LocalVar(ins)
ham_detuning = Empty
for i in range(n):
    ham_detuning += -d * noper[i]
ins.set_ham(ham_detuning)

# for i in range(n):
L = SignalLine(Rydberg)
ins = Instruction(L, "native", f"Rabi")
o = LocalVar(ins)
p = LocalVar(ins)
ham_rabi = Empty
for i in range(n):
    ham_rabi += o / 2 * (Expression.cos(p) * q[i].X - Expression.sin(p) * q[i].Y)
ins.set_ham(ham_rabi)

# ins.set_ham(o / 2 * (Expression.cos(p) * q[i].X - Expression.sin(p) * q[i].Y))
