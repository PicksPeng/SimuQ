import numpy as np

from simuq.environment import Qubit
from simuq.expression import Expression
from simuq.hamiltonian import hlist_sum
from simuq.qmachine import QMachine


C_6 = 862690 * 2.0 * np.pi


def generate_qmachine(n=3, inits=None):
    rydberg = QMachine()

    q = [Qubit(rydberg) for i in range(n)]

    l = C_6 ** (1.0 / 6)
    if inits is None:
        x = [0] + [rydberg.add_global_variable(init_value=l * i) for i in range(1, n)]
    else:
        x = [0] + [rydberg.add_global_variable(init_value=inits[i]) for i in range(1, n)]

    noper = [(q[i].I - q[i].Z) / 2 for i in range(n)]

    hlist = []
    for i in range(n):
        for j in range(i):
            hlist.append((C_6 / (x[i] - x[j]) ** 6) * noper[i] * noper[j])
    sys_h = hlist_sum(hlist)
    rydberg.set_sys_ham(sys_h)

    L = rydberg.add_signal_line()
    ins = L.add_instruction("native", "Detuning")
    d = ins.add_local_variable()
    ham_detuning = 0
    for i in range(n):
        ham_detuning += -d * noper[i]
    ins.set_ham(ham_detuning)

    L = rydberg.add_signal_line()
    ins = L.add_instruction("native")
    o = ins.add_local_variable()
    p = ins.add_local_variable()
    ham_rabi = 0
    for i in range(n):
        ham_rabi += o / 2 * (Expression.cos(p) * q[i].X - Expression.sin(p) * q[i].Y)
    ins.set_ham(ham_rabi)

    return rydberg
