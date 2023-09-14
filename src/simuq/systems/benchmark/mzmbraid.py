# A tri-junction Majorana zero mode braiding experiment
# The following is decribed in [Stenger et al., Phys. Rev. Research 3 033171, 2021]

from numpy import linspace

from simuq.environment import Fermion
from simuq.qsystem import QSystem


def GenQS(tau=3.3, alpha=3, Jmax=1, D=3, segs=6):
    qs = QSystem()
    f = [Fermion(qs) for i in range(3)]

    gamma_x = [f[i].a + f[i].c for i in range(3)]
    gamma_y = [-1j * (f[i].a - f[i].c) for i in range(3)]

    def model(alpha, J01, J12, J20):
        J = [[0, J01, -J20], [-J01, 0, J12], [J20, -J12, 0]]
        h = 0
        for i in range(3):
            h += 1j * alpha * gamma_x[i] * gamma_y[i]
        for i in range(3):
            for j in range(3):
                h += 0.5j * J[i][j] * gamma_x[i] * gamma_x[j]
        return h

    def branch(J01_init, J01_end, J12_init, J12_end, J20_init, J20_end):
        def ret(t):
            J01 = J01_init + (J01_end - J01_init) * (t / tau)
            J12 = J12_init + (J12_end - J12_init) * (t / tau)
            J20 = J20_init + (J20_end - J20_init) * (t / tau)
            return model(alpha, J01, J12, J20)

        return ret

    params = [[Jmax, 0, 0, Jmax, 0, 0], [0, 0, Jmax, 0, 0, Jmax], [0, Jmax, 0, 0, Jmax, 0]]
    params = params + params

    for i in range(segs):
        qs.add_td_evolution(branch(*params[i]), linspace(0, tau, D + 1))

    return qs
