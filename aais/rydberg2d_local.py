from simuq.environment import qubit
from simuq.qmachine import *
from simuq.expression import Expression
import numpy as np

def GenMach(n = 3) :
    Rydberg = QMachine()

    C6 = 862690 * 2. * np.pi

    q = [qubit(Rydberg) for i in range(n)]

    l = (C6 / 4) ** (1. / 6) / (2 - 2 * np.cos(2 * np.pi / n)) ** 0.5
    x = [(0, 0.)] + [(GlobalVar(Rydberg, init_value = l * (np.cos(i * 2 * np.pi / n) - 1)), GlobalVar(Rydberg, init_value = l * np.sin(i * 2 * np.pi / n))) for i in range(1, n)]
    noper = [(q[i].I - q[i].Z) / 2 for i in range(n)]

    sys_h = 0
    for i in range(n) :
        for j in range(i) :
            dsqr = (x[i][0] - x[j][0])**2 + (x[i][1] - x[j][1])**2
            sys_h += (C6 / (dsqr ** 3)) * noper[i] * noper[j]
    Rydberg.set_sys_ham(sys_h)

    L = SignalLine(Rydberg)
    ins = Instruction(L, "native", f"Detuning")
    d = LocalVar(ins)
    ham_detuning = 0
    for i in range(n):
        ham_detuning += -d * noper[i]
    ins.set_ham(ham_detuning)

    L = SignalLine(Rydberg)
    ins = Instruction(L, "native", f"Rabi")
    o = LocalVar(ins)
    p = LocalVar(ins)
    ham_rabi = 0
    for i in range(n):
        ham_rabi += o / 2 * (Expression.cos(p) * q[i].X - Expression.sin(p) * q[i].Y)
    ins.set_ham(ham_rabi)

    return Rydberg
