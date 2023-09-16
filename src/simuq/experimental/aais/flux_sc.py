import numpy as np

from simuq.environment import Qubit
from simuq.expression import Expression
from simuq.qmachine import *

FluxSCMach = QMachine()
mach = FluxSCMach
L = 5
ql = [Qubit(mach) for i in range(L)]
Js = [
    2 * np.pi * 0.01084,
    2 * np.pi * 0.00028,
]  # We assume they are transverse here, which can be specified in detail.
l = len(Js)


for i in range(L):
    Line = SignalLine(mach)

    ins1 = Instruction(Line, "native", "L{}_X_Y".format(i))
    amp = LocalVar(ins1)
    phase = LocalVar(ins1)
    ins1.set_ham(amp * (Expression.cos(phase) * ql[i].X + Expression.sin(phase) * ql[i].Y))

    ins2 = Instruction(Line, "derived", "L{}_Z".format(i))
    amp = LocalVar(ins2)
    ins2.set_ham(amp * ql[i].Z)

Line = SignalLine(mach)

ins = Instruction(Line, "native", "L_int")
amp = LocalVar(ins)
hint = 0
for j in range(1, l + 1):
    for i in range(L - j):
        hint += Js[j - 1] / 2 * (ql[i].X * ql[i + j].X + ql[i].Y * ql[i + j].Y)
ins.set_ham(amp * hint)
