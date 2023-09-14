import numpy as np

from simuq.environment import Qubit
from simuq.expression import Expression
from simuq.hamiltonian import TIHamiltonian
from simuq.qmachine import QMachine


def ham_sum(hlist):
    n = len(hlist)
    if n == 0:
        return 0
    sites_type = hlist[0].sites_type
    for i in range(n):
        if hlist[i].sites_type != sites_type:
            raise Exception("Site types do not match")
    ham = []
    for i in range(n):
        ham += hlist[i].ham

    return TIHamiltonian(sites_type, ham)


C_6 = 862690 * 2.0 * np.pi


def generate_qmachine(n=3, inits=None):
    rydberg = QMachine()

    q = [Qubit(rydberg) for i in range(n)]

    l = (C_6 / 4) ** (1.0 / 6) / (2 - 2 * np.cos(2 * np.pi / n)) ** 0.5

    if inits is None:
        x = [(0.0, 0.0)] + [
            (
                rydberg.add_global_variable(init_value=l * (np.cos(i * 2 * np.pi / n) - 1)),
                rydberg.add_global_variable(init_value=l * np.sin(i * 2 * np.pi / n)),
            )
            for i in range(1, n)
        ]
    else:
        x = [(0.0, 0.0)] + [
            (
                rydberg.add_global_variable(init_value=inits[i][0]),
                rydberg.add_global_variable(init_value=inits[i][1]),
            )
            for i in range(1, n)
        ]
    noper = [(q[i].I - q[i].Z) / 2 for i in range(n)]

    hlist = []
    for i in range(n):
        for j in range(i):
            dsqr = (x[i][0] - x[j][0]) ** 2 + (x[i][1] - x[j][1]) ** 2
            hlist.append((C_6 / (dsqr**3)) * noper[i] * noper[j])
    sys_h = ham_sum(hlist)

    rydberg.set_sys_ham(sys_h)

    for i in range(n):
        L = rydberg.add_signal_line()
        ins = L.add_instruction("native", f"Detuning of site {i}")
        d = ins.add_local_variable()
        ins.set_ham(-d * noper[i])

    for i in range(n):
        L = rydberg.add_signal_line()
        ins = L.add_instruction("native", f"Rabi of site {i}")
        o = ins.add_local_variable()
        p = ins.add_local_variable()
        ins.set_ham(o / 2 * (Expression.cos(p) * q[i].X - Expression.sin(p) * q[i].Y))

    return rydberg
