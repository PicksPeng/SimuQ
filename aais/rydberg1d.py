from simuq.environment import qubit, fock
from simuq.qmachine import *
from simuq.expression import Expression
from simuq.hamiltonian import Empty
import numpy as np

def GenMach(n = 3) :
    Rydberg = QMachine()

    C6 = 862690 * 2 * np.pi

    q = [qubit(Rydberg) for i in range(n)]

    l = C6 ** (1. / 6)
    x = [0] + [GlobalVar(Rydberg, init_value = l * i) for i in range(1, n)]
    noper = [(q[i].I - q[i].Z) / 2 for i in range(n)]

    sys_h = Empty
    for i in range(n) :
        for j in range(i) :
            sys_h += (C6 / (x[i] - x[j])**6) * noper[i] * noper[j]
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

    return Rydberg
