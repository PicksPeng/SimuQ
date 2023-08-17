from simuq.aais.qmachine_factory import QMachineFactory
from simuq.qmachine import *
from simuq.expression import Expression
from simuq.hamiltonian import TIHamiltonian
import numpy as np

def ham_sum(hlist) :
    n = len(hlist)
    if n == 0 :
        return 0
    sites_type = hlist[0].sites_type
    for i in range(n) :
        if hlist[i].sites_type != sites_type :
            raise Exception("Site types do not match")
    ham = []
    for i in range(n) :
        ham += hlist[i].ham

    return TIHamiltonian(sites_type, ham)

class Rydberg1DQMachineFactory(QMachineFactory):
    @staticmethod
    def generate_qmachine(n=3, inits=None, *args, **kwargs):
        rydberg = QMachine()

        C6 = 862690 * 2. * np.pi

        q = [qubit(rydberg) for i in range(n)]

        l = C6 ** (1. / 6)
        if inits == None :
            x = [0] + [GlobalVar(rydberg, init_value = l * i) for i in range(1, n)]
        else :
            x = [0] + [GlobalVar(rydberg, init_value = inits[i]) for i in range(1, n)]

        noper = [(q[i].I - q[i].Z) / 2 for i in range(n)]

        hlist = []
        for i in range(n) :
            for j in range(i) :
                hlist.append((C6 / (x[i] - x[j])**6) * noper[i] * noper[j])
        sys_h = ham_sum(hlist)
        rydberg.set_sys_ham(sys_h)

        for i in range(n) :
            L = SignalLine(rydberg)
            ins = Instruction(L, 'native', f'Detuning of site {i}')
            d = LocalVar(ins)
            ins.set_ham(- d * noper[i])

        for i in range(n) :
            L = SignalLine(rydberg)
            ins = Instruction(L, 'native')
            o = LocalVar(ins)
            p = LocalVar(ins)
            ins.set_ham(o / 2 * (Expression.cos(p) * q[i].X - Expression.sin(p) * q[i].Y))

        return rydberg
